#include "BlockCube.h"
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

BlockCube::BlockCube()
{
}

BlockCube::BlockCube(const BlockCube&)
{
}

BlockCube::~BlockCube()
{
}

BlockCube&
BlockCube::operator=(const BlockCube&)
{
        return *this;
}

void BlockCube::QueryCube(std::vector<int> query, std::string queries_ops, std::map<std::vector<int>, int>& query_cache, int my_rank, int num_dims, int num_tuples, std::string output_folder, int num_procs, int tuple_partition_size, bool silent){

	//Verifica no cache de consultas se a consulta já foi repetida para evitar retrabalho - útil nas inquires
	if(query_cache.count(query) == 0){

	//Cada processo MPI tem um diretório próprio para armazenar dados do cubo
	std::string process_directory = output_folder + "/proc" + std::to_string(my_rank);

	//Listas de BIDs associados a valores de atributo da consulta
	//Cada posição do vector armazena uma lista de BIDs
	std::vector<std::vector<int>> lists_of_bids;

	//Listas de TIDs fazem parte da resposta da consulta
	//Cada posição do vector armazena uma lista de TIDs
	std::vector<std::vector<int>> lists_of_tids;

	//Quantidade de BIDs criados para o cubo, na partição de tuplas recebida
	int num_bids = (tuple_partition_size % tbloc == 0) ? tuple_partition_size / tbloc : tuple_partition_size / tbloc + 1;

	//Lista de BIDs, após ter feito a interseção de todas com base nas points
	std::vector<int> bids_intersection;

	//Lista final com a interseção de todos os TIDs que respondem à consulta
	std::vector<int> tids_intersection;

	//No caso do inquire irá armazenar listas de valores de atributo
	//Essas listas serão usadas para gerar as consultas pontuais derivadas do inquire
	std::vector<std::vector<int>> lists_of_attribs(num_dims);

	//Lista de dimensões sobre as quais incide um operador inquire
	std::vector<int> inquires;

	//Lista de dimensões sobre as quais incide um operador point
	std::vector<int> points;

	//Lista de dimensões sobre as quais incide um operador agregação
	std::vector<int> aggregations;

	//No momento a consulta retorna apenas COUNT, então aqui é armazenado o resultado local
	int local_count = 0;

	//Essa variável só deve ser usada pelo processo de menor rank para guardar o resultado final
	int global_count = 0;

	//O primeiro passo para responder uma consulta com o algoritmo bCubing é resolver operandos do tipo point
	for(std::vector<T>::size_type operand = 0; operand != query.size(); operand++) {

		//Cada operando está associado a uma das dimensões
		//O primeiro operando da consulta é da primeira dimensão, o segundo da segunda, etc.
		int dim = operand;

		//Se o operando não for do tipo agregação '*' ou inquire '?', então é do tipo point, um valor de atributo
		if (query[operand] != -1 && query[operand] != -2)
		{
			//Extrai a lista de BIDs, aqui ainda no formato de um conjunto sem duplicações
			std::set<int> bids_set = bblocRAM[dim][query[operand]];

			//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
			std::vector<int> bids_vector(bids_set.begin(), bids_set.end());

			//Adiciona às listas de BIDs a lista de BIDs do valor de atributo buscando na bblocRAM
			lists_of_bids.push_back(bids_vector);

			//Atualiza a lista de dimensões afetada por points
			points.push_back(dim);
		}

		//Se o operando for do tipo inquire '?', identificamos a dimensão afetada
		if (query[operand] == -2){
			//Atualiza a lista de dimensões afetada por inquires
			inquires.push_back(dim);
		}

		//Se o operando for do tipo agregação '*', identificamos a dimensão afetada
		if (query[operand] == -1){
			//Atualiza a lista de dimensões afetada por inquires
			aggregations.push_back(dim);
		}

	}

	//Se não tiver nenhuma point, devem ser todas agregações ou podem ter inquires!
	if(points.size() == 0){

        //Se só há inquires ou agregações na consulta, então todos os BIDs devem ser usados
        for(int bid = 0; bid < num_bids; bid++){
        	bids_intersection.push_back(bid);
        }

	} else {//Caso contrário há points na consulta

		//Faz a interseção das listas de BIDS encontradas
		bids_intersection = IntersectMultipleVectors(lists_of_bids);
	}

	//Se não tiver nenhuma consulta inquire, deve ter apenas points
	if(inquires.empty()){

		//Só faz sentido tentar responder a consulta se a interseção não for fazia
		if(!bids_intersection.empty()){

			//Se todas as operações da consulta forem agregações, temos o caso: * * * *
			if(aggregations.size() == static_cast<std::vector<int>::size_type>(num_dims)){

				//Soma de tuplas até o começo da partição associada a esse processo
				int sum_of_partitions = std::accumulate(tuple_partition_listings.begin(), tuple_partition_listings.begin() + my_rank, 0);

				//Adiciona TODOS os valores de TIDs da partição
				for(int tid = 1; tid <= tuple_partition_size; tid++){
					tids_intersection.push_back(tid + sum_of_partitions);
				}

				//Atualiza o valor de COUNT localmente para a consulta
				local_count = tids_intersection.size();

				//Adicionalmente, verificamos se tem alguma outra operação de agregação diferente de Frequência (Fr)
				if(!queries_ops.empty()){

					for(auto & dim_number : aggregations){

						//Se tiver precisamos buscar as medidas de disco da mesma forma
						for(auto & bid : bids_intersection){

							//Esse é o diretório do bloco de BID associado aos dados que desejamos, já existente em disco
							std::string bloc_directory = process_directory + "/dim" + std::to_string(dim_number) + "/bloc" + std::to_string(bid);

							//LÊ AS MEDIDAS DO BLOCO EM DISCO

							//Nome do arquivo onde as medidas do BID estão salvas - diretório do processo, num diretório específico da dimensão
							std::string bid_meas_filename = bloc_directory + "/meas/" + std::to_string(bid) + ".bin";

							//Indica que o arquivo de entrada será um stream de dados binários
							std::ifstream ifm(bid_meas_filename.c_str(), std::ifstream::binary);

							//Serialização é feita pela biblioteca Boost
							boost::archive::binary_iarchive im(ifm, boost::archive::no_header);

							//Carrega os dados de medidas do BID que estavam em disco
							im & bmeas[bid];

						}
					}
				}

			} else { //Necessariamente deve ter alguma point nessa consulta

				//Se não tiver outras operações, então deve fazer apenas a verificação de frequência (Fr)
				//Nesse caso não precisa acessar disco pois já fizemos o cache em memória principal, porém só funciona ]
				//se tiver apenas uma point
				if(queries_ops.empty() && points.size() == 1){

					//Para cada uma das dimensões associadas à consultas point
					for(auto & dim_number : points){

						//Para cada um dos BIDs da interseção
						for(auto & bid : bids_intersection){

							//Acessamos o cache de frequência em memória, somando a frequência em cada BID
							local_count += bblocRAMFreq[dim_number][bid][query[dim_number]];

						}

					}

				} else { //Tem outras medidas a serem computadas, então deve buscar TIDs de disco

					//Agora sabemos quantos BIDs teremos que recuperar do disco para responder a consulta
					for (auto & dimension : bbloc) {
						//Redimensiona a estrutura em memória da bCubingBloc para acomodar os blocos em cada dimensão
						//Como usamos o índice para determinar o BID, e a interseção de BIDs está ordenada, o último valor deve ser o maior índice
						//Assim alocamos memória suficiente para todos os índices até o maior, muitos índicies ficaram vazios
						//
						// E.g. bids_intersection = [1, 3]
						//
						// Então irá alocar do índice 0 até o 3, inicialmente vazio, na bbloc [[ ],[ ],[ ],[ ]]
						// Depois leremos de disco para preencher conforme o necessário
						//
						dimension.resize(bids_intersection.back() + 1);
					}

					//Para cada uma das dimensões associadas à consultas point
					for(auto & dim_number : points){

						//Armazena os TIDs associados ao valor de atributo de uma das dimensões da consulta
						std::vector<int> tids;

						//Resolveremos primeiro elas, buscando os BIDs associados em disco
						for(auto & bid : bids_intersection){

							//LÊ O BLOCO DO DISCO

							//Esse é o diretório do bloco de BID associado aos dados que desejamos, já existente em disco
		                    std::string bloc_directory = process_directory + "/dim" + std::to_string(dim_number) + "/bloc" + std::to_string(bid);

							//Nome do arquivo onde o BID está salvo - diretório do processo, num diretório específico da dimensão
							std::string bid_filename = bloc_directory + "/tids/" + std::to_string(bid) + ".bin";

							//Indica que o arquivo de entrada será um stream de dados binários
							std::ifstream ifb(bid_filename.c_str(), std::ifstream::binary);

							//Serialização é feita pela biblioteca Boost
							boost::archive::binary_iarchive ib(ifb, boost::archive::no_header);

							//Carrega os dados do BID que estava em disco
							ib & bbloc[dim_number][bid];

							//Acessa uma única vez a estrutura bCubingBloc para obter os TIDS para esse BID
							std::vector<int> tids_bbloc = bbloc[dim_number][bid][query[dim_number]];

							//Adiciona às listas de TID os TIDs da dimensão associada à point, para esse BID, com base no valor da consulta
							tids.insert(std::end(tids), std::begin(tids_bbloc), std::end(tids_bbloc));

		                    //LÊ AS MEDIDAS DO BLOCO EM DISCO

		                    //Nome do arquivo onde as medidas do BID estão salvas - diretório do processo, num diretório específico da dimensão
		                    std::string bid_meas_filename = bloc_directory + "/meas/" + std::to_string(bid) + ".bin";

		                    //Indica que o arquivo de entrada será um stream de dados binários
		                    std::ifstream ifm(bid_meas_filename.c_str(), std::ifstream::binary);

		                    //Serialização é feita pela biblioteca Boost
		                    boost::archive::binary_iarchive im(ifm, boost::archive::no_header);

		                    //Carrega os dados de medidas do BID que estavam em disco
		                    im & bmeas[bid];

						}

						//Vai armazenar uma quantidade de listas igual à quantidade de consultas do tipo point
						lists_of_tids.push_back(tids);

					}

					//Faz a interseção das listas de TIDS encontradas
					tids_intersection = IntersectMultipleVectors(lists_of_tids);

					//Atualiza o valor de COUNT localmente para a consulta
					local_count = tids_intersection.size();

				}
			}

		}

		//Somente será executado quando todos processos chegarem nesse ponto
		//Combina as respostas de todos, e todos sabem qual a contagem global de tuplas
		MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//Se essa consulta teve algum resultado, vamos guardá-la em cache para evitar a reconsulta
		if(global_count > 0){
			query_cache[query] = global_count;
		}

		if(my_rank == 0){
			if(global_count > 0){
				if(!silent){
					for(auto & operand : query) {
						(operand == -1) ? std::cout << '*' << ' ' : std::cout << operand << ' ';
					}
					std::cout << ": " << global_count << ' ';
				}
			}
		}

		//Agora já fizemos a contagem, que é o padrão, mas podem ter mais operações de agregação
		if(!queries_ops.empty()){

	        //Armazena um dos operadores fornecidos
	        std::string arg;

	        //Guarda um identificador da medida associada, a primeira medida é M0
	        int measure_id = 0;

	        //A lista de operações, em formato de string
	        std::stringstream ops(queries_ops);

	        //Passa por cada argumento da lista de operações, separado por espaços
	        while (ops >> arg) {

	        	//Operação de agração do tipo SOMA
	        	if(arg == "Sm" || arg == "soma"){

	        		//A soma local a princípio é zero
	        		float local_sum = 0;

	        		//Essa variável só deve ser usada pelo processo de menor rank para guardar o resultado final
	        		float global_sum = 0;

	        		//Caso a consulta tenha resultado em alguns tids
	        		if(local_count > 0){

	        			//Soma localmente os valores com base nos TIDs da resposta
						for(auto & tid : tids_intersection){

							//No caso do bCubing as medidas são indexadas pelo BID
							int bid;

							//Soma de tuplas até o começo da partição associada a esse processo
							int sum_of_partitions = std::accumulate(tuple_partition_listings.begin(), tuple_partition_listings.begin() + my_rank, 0);

							//Aqui fazemos o cálculo para saber qual o BID o TID que queremos está
							if(my_rank > 0){
								//Se não for o processo de menor rank, precisamos determinar em que ponto o TID está na partição do processo
								//Por exemplo, se a partição do processo começa no TID 4, então o TID 4 para ele é equivalente ao TID 1 do primeiro processo
								//Ambos os TIDs estariam no primeiro bloco de cada processo
								//Ex: se TID 4 então 4 % ( 3 + 1 ) + 1 = 4 % (4) + 1 = 0 + 1 = 1
								//Ou seja, os processos anteriores foram até o TID 3 e nesse processo a partição começa no TID 4 que para ele é o primeiro TID
								//Ex: tbloc = 2 e tid = 1 então bid = 1/2 + 1%2 - 1 = 0 + 1 - 1 = 0, está no BID 0
								bid = (tid % (sum_of_partitions + 1) + 1)/tbloc + ((tid % (sum_of_partitions + 1) + 1) % tbloc != 0) - 1;
							} else {
								//Se for o processo de menor rank, então ele pegou a primeira partição
								//Nesse caso dá pra saber o BID dividindo o valor dele pelo tamanho do bloco e arrendondando pra cima
								//Ex: tbloc = 2 e tid = 3 então bid = 3/2 + 3%2 - 1 = 1 + 1 - 1 = 1, está no BID 1
								bid = tid/tbloc + (tid % tbloc != 0) - 1;
							}

							//Acessa o índice em memória para obter o valor de medida e faz a SOMA
							local_sum += bmeas[bid][tid][measure_id];
						}

	        		}

	        		//Somente será executado quando todos processos chegarem nesse ponto
	        		//Combina as respostas de todos num só processo
	        		MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	        		//Apresenta o valor agregado de SOMA das medidas
	        		if(my_rank == 0){
	        			if(global_count > 0){
	        				if(!silent){
		        				std::cout << "Sm(M" << measure_id << ") = " << global_sum << ' ';
	        				}
	        			}
	        		}

	        	}

	        	//Operação de agregação do tipo MEDIANA
	        	if(arg == "Md" || arg == "mediana"){

	        		//Como é uma operação holística deve apenas extrair as medidas inicialmente
	        		std::vector<float> local_measures;

	        		//Essa variável só deve ser usada pelo processo de menor rank para guardar o resultado final
	        		std::vector<float> global_measures(global_count);

	        		//Lista de quantidades de resultados obtidos por cada processo para a consulta
	        		std::vector<int> recvcounts(num_procs);

	        		//Lista de deslocamentos necessários para organizar todos os resultados numa única lista global
	        		std::vector<int> displs(num_procs);

	        		//Caso a consulta tenha resultado em alguns tids
	        		if(local_count > 0){

	        			//Soma localmente os valores com base nos TIDs da resposta
						for(auto & tid : tids_intersection){

							//No caso do bCubing as medidas são indexadas pelo BID
							int bid;

							//Soma de tuplas até o começo da partição associada a esse processo
							int sum_of_partitions = std::accumulate(tuple_partition_listings.begin(), tuple_partition_listings.begin() + my_rank, 0);

							//Aqui fazemos o cálculo para saber qual o BID o TID que queremos está
							if(my_rank > 0){
								//Se não for o processo de menor rank, precisamos determinar em que ponto o TID está na partição do processo
								//Por exemplo, se a partição do processo começa no TID 4, então o TID 4 para ele é equivalente ao TID 1 do primeiro processo
								//Ambos os TIDs estariam no primeiro bloco de cada processo
								//Ex: se TID 4 então 4 % ( 3 + 1 ) + 1 = 4 % (4) + 1 = 0 + 1 = 1
								//Ou seja, os processos anteriores foram até o TID 3 e nesse processo a partição começa no TID 4 que para ele é o primeiro TID
								//Ex: tbloc = 2 e tid = 1 então bid = 1/2 + 1%2 - 1 = 0 + 1 - 1 = 0, está no BID 0
								bid = (tid % (sum_of_partitions + 1) + 1)/tbloc + ((tid % (sum_of_partitions + 1) + 1) % tbloc != 0) - 1;
							} else {
								//Se for o processo de menor rank, então ele pegou a primeira partição
								//Nesse caso dá pra saber o BID dividindo o valor dele pelo tamanho do bloco e arrendondando pra cima
								//Ex: tbloc = 2 e tid = 3 então bid = 3/2 + 3%2 - 1 = 1 + 1 - 1 = 1, está no BID 1
								bid = tid/tbloc + (tid % tbloc != 0) - 1;
							}

							//Obtém a lista de medidas localmente
							local_measures.push_back(bmeas[bid][tid][measure_id]);
						}

	        		}

	        		//Compartilha com os demais processos quantos resultados obteve localmente
	                MPI_Allgather(&local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

	                //Com base na quantidade de resultados dos demais, o processo consegue saber o deslocamento necessário
	                //para armazenar os seus próprios resultados na lista que o  processo principal irá gerar com todos os resultados
	                //Ex: processo 0 tem 3 resultados, se eu sou o processo 1 meus resultados na lista final deve começar a partir do 3
	                int displacement = std::accumulate(recvcounts.begin(), recvcounts.begin() + my_rank, 0);

	                //Compartilha com os demais processos os deslocamentos
	                MPI_Allgather(&displacement, 1, MPI_INT, displs.data(), 1, MPI_INT, MPI_COMM_WORLD);

	                //O processo principal de rank 0 coleta as medidas dos demais considerando os deslocamentos necessários
	                //Os resultados são salvos em uma lista global de medidas com as medidas obtidas por todos os processos
	        		MPI_Gatherv(local_measures.data(), local_measures.size(), MPI_FLOAT,
	        				global_measures.data(), recvcounts.data(), displs.data(), MPI_FLOAT,
	        		    0, MPI_COMM_WORLD);

	        		//Apresenta o valor agregado de MEDIANA das medidas
	        		if(my_rank == 0){
	        			if(global_count > 0){

	        				//Usamos a estratégia de ordenar apenas o suficiente a lista de medidas
	        			    size_t n = global_measures.size() / 2;

	        			    //Para obter a mediana é preciso que a lista esteja ordenada, nesse caso apenas o suficiente
	        			    std::nth_element(global_measures.begin(), global_measures.begin()+n, global_measures.end());

	        			    if(!silent){
		        			    //Obtemos o valor mediano, após a ordenação
		        				std::cout << "Md(M" << measure_id << ") = " << global_measures[n] << ' ';
	        			    }
	        			}
	        		}

	        	}

	        	//O próximo operador será associado à próxima medida
	        	measure_id++;

	        }
		}

		//Final de linha, apenas para organizar as respostas
		if(my_rank == 0){
			if(global_count > 0){
				if(!silent){
					//Oficialmente terminamos de responder a consulta!
					std::cout << std::endl;
				}
			}
		}

	} else { //Consulta tem pelo menos um operador inquire

		//Para cada uma das dimensões associadas à consultas point
		for(auto & dim_number : points){
			//Insere o valor de atributo como uma nova lista nas listas de atributos
			lists_of_attribs[dim_number] = { query[dim_number] };
		}

		//Para cada uma das dimensões associadas à agregações
		for(auto & dim_number : aggregations){
			//Insere o valor de atributo como uma nova lista nas listas de atributos
			lists_of_attribs[dim_number] = { query[dim_number] };
		}

		//Só faz sentido tentar responder a consulta se a interseção não for fazia
		if(!bids_intersection.empty()){

			//Agora sabemos quantos BIDs teremos que recuperar do disco para responder a consulta
			for (auto & dimension : bbloc) {
				//Redimensiona a estrutura em memória da bCubingBloc para acomodar os blocos em cada dimensão
				//Como usamos o índice para determinar o BID, e a interseção de BIDs está ordenada, o último valor deve ser o maior índice
				//Assim alocamos memória suficiente para todos os índices até o maior, muitos índicies ficaram vazios
				//
				// E.g. bids_intersection = [1, 3]
				//
				// Então irá alocar do índice 0 até o 3, inicialmente vazio, na bbloc [[ ],[ ],[ ],[ ]]
				// Depois leremos de disco para preencher conforme o necessário
				//
				dimension.resize(bids_intersection.back() + 1);
			}

			//Para cada uma das dimensões associadas à consultas inquire
			for(auto & dim_number : inquires){

				//Irá guardar os valores de atributo possíveis na dimensão
				std::set<int> attribs_set;

				//Um valor extra necessário é o de agregação '*' pois não está nos dados e é usado em consulta válida de inquire
				//E.g: uma dimensão com atributos [1,2] na verdade tem os atributos [-1,1,2], onde -1 indica TODOS
				attribs_set.insert(-1);

				//Resolveremos primeiro elas, buscando os BIDs associados em disco
				for(auto & bid : bids_intersection){

					//LÊ O BLOCO DO DISCO

					//Esse é o diretório do bloco de BID associado aos dados que desejamos, já existente em disco
                    std::string bloc_directory = process_directory + "/dim" + std::to_string(dim_number) + "/bloc" + std::to_string(bid);

					//Nome do arquivo onde o BID está salvo - diretório do processo, num diretório específico da dimensão
					std::string bid_filename = bloc_directory + "/tids/" + std::to_string(bid) + ".bin";

					//Indica que o arquivo de entrada será um stream de dados binários
					std::ifstream ifb(bid_filename.c_str(), std::ifstream::binary);

					//Serialização é feita pela biblioteca Boost
					boost::archive::binary_iarchive ib(ifb, boost::archive::no_header);

					//Carrega os dados do BID que estava em disco
					ib & bbloc[dim_number][bid];

					//Para cada par associando um valor de atributo a uma lista de TIDs
					for (auto& attrib_pair_tids: bbloc[dim_number][bid]) {
						//Extraia apenas o valor de atributo, o primeiro elemento do pair, e salve na lista de atributos
						attribs_set.insert(attrib_pair_tids.first);
					}

                    //LÊ AS MEDIDAS DO BLOCO EM DISCO

                    //Nome do arquivo onde as medidas do BID estão salvas - diretório do processo, num diretório específico da dimensão
                    std::string bid_meas_filename = bloc_directory + "/meas/" + std::to_string(bid) + ".bin";

                    //Indica que o arquivo de entrada será um stream de dados binários
                    std::ifstream ifm(bid_meas_filename.c_str(), std::ifstream::binary);

                    //Serialização é feita pela biblioteca Boost
                    boost::archive::binary_iarchive im(ifm, boost::archive::no_header);

                    //Carrega os dados de medidas do BID que estavam em disco
                    im & bmeas[bid];

				}

				//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
				std::vector<int> attribs_vector(attribs_set.begin(), attribs_set.end());

				//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
				lists_of_attribs[dim_number] = attribs_vector;

			}


		}

		//Nesse ponto cada processo deve ter gerado listas de listas de atributos
		//Dimensões onde há inquire terão listas maiores, algo do tipo para a consulta ? 1 2 = [[1,2,3,4,5], 1, 2]
		//Então o que eu quero agora é: cada processo a partir do 0, para cada inquire, vai enviar os valores de atributos que tem
		for(auto & dim_number : inquires){

			//Esse set irá guardar os atributos, sem repetições, entre os obtidos por todos os processos
			std::set<int> proc_attribs_inquire_unique;

			//Para cada uma das dimensões associadas à consultas inquire
			for(int proc = 1; proc < num_procs; proc++){

				//Número de atributos que serão enviados
				int num_attribs;

				//Irá armazenar os valores de atributo de cada processo, para cada inquire
				std::vector<int> proc_attribs_inquire;

				//Se for o rank 0 irá receber os atributos
				if(my_rank == 0){

					//Primeiro precisa receber um inteiro com o tamanho da lista que será enviada
				    MPI_Recv(&num_attribs, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				    //Garanta que o vector tem espaço suficiente
				    proc_attribs_inquire.resize(num_attribs);

					//Recebe a lista de atributos
				    MPI_Recv(proc_attribs_inquire.data(), num_attribs, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				    //Coloque todos os atributos do vector no set, sem repetições
				    for(auto &attr : proc_attribs_inquire){
				    	proc_attribs_inquire_unique.insert(attr);
				    }

				} else if (my_rank == proc){ //Demais irão enviar seus atributos (o 0 envia para ele mesmo)

					//Calcula o tamanho da lista que será enviada, com os atributos dessa dimensão com inquire
					num_attribs = lists_of_attribs[dim_number].size();

					//Envia o tamanho da lista de atributos da dimensão com inquire
				    MPI_Send(&num_attribs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

					//Envia a lista de atributos
				    MPI_Send(lists_of_attribs[dim_number].data(), num_attribs, MPI_INT, 0, 0, MPI_COMM_WORLD);

				}
			}

			//Aqui o processo de menor rank deve ter uma lista com todos os atributos!
			//Só faltam os atributos dele próprio, entao fazemos isto
		    //Coloque todos os atributos do vector no set, sem repetições
		    for(auto &attr : lists_of_attribs[dim_number]){
		    	proc_attribs_inquire_unique.insert(attr);
		    }

			//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
			std::vector<int> attribs_vector(proc_attribs_inquire_unique.begin(), proc_attribs_inquire_unique.end());

			//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
			lists_of_attribs[dim_number] = attribs_vector;

		}

		//Agora o processo de menor rank tem uma lista de todos os atributos de cada inquire
		//Fazemos o processo inverso, atualizando os demais processo com essas listas
		for(auto & dim_number : inquires){

			//Para cada uma das dimensões associadas à consultas inquire
			for(int proc = 1; proc < num_procs; proc++){

					//Número de atributos que serão enviados
					int num_attribs;

					//Irá armazenar os valores de atributo de cada processo, para cada inquire
					std::vector<int> proc_attribs_inquire;

					//Se não for o rank 0 irá receber os atributos
					if(my_rank == proc){

						//Primeiro precisa receber um inteiro com o tamanho da lista que será enviada
						MPI_Recv(&num_attribs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						//Garanta que o vector tem espaço suficiente
						proc_attribs_inquire.resize(num_attribs);

						//Recebe a lista de atributos
						MPI_Recv(proc_attribs_inquire.data(), num_attribs, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
						lists_of_attribs[dim_number] = proc_attribs_inquire;


					} else if (my_rank == 0){ //O 0 envia seus atributos para os demais

						//Calcula o tamanho da lista que será enviada, com os atributos dessa dimensão com inquire
						num_attribs = lists_of_attribs[dim_number].size();

						//Envia o tamanho da lista de atributos da dimensão com inquire
						MPI_Send(&num_attribs, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);

						//Envia a lista de atributos
						MPI_Send(lists_of_attribs[dim_number].data(), num_attribs, MPI_INT, proc, 0, MPI_COMM_WORLD);

					}
			}

		}

		//Agora serão geradas as combinações de consultas associadas ao inquire
		//https://www.geeksforgeeks.org/combinations-from-n-arrays-picking-one-element-from-each-array/

		//Conta a quantidade de listas de atributos
		int number_of_lists = lists_of_attribs.size();

		//Guarda posições de índices de cada uma das listas
		//Essas posições indicam qual o índice da lista está sendo usado para a combinação
		int* indices = new int[number_of_lists];

		//Inicialmente as combinações são geradas partindo do primeiro índice das listas
		for (int i = 0; i < number_of_lists; i++){
			indices[i] = 0;
		}

		//Só sai desse laço com um "break"
		while (1) {

			//Armazena a consulta gerada pelo processo nessa iteração
			std::vector<int> query(num_dims);

			//Armazena a consulta, dentre aquelas gerados por todos os processos nessa iteração, sendo executada nesse momento
			std::vector<int> query_running(num_dims);

			//Gera a consulta com base na próxima combinaçao de atributos
			for (int i = 0; i < number_of_lists; i++){
				query[i] = lists_of_attribs[i][indices[i]];
			}

			//Todos os processos executam a query
			QueryCube(query, queries_ops, query_cache, my_rank, num_dims, num_tuples, output_folder, num_procs, tuple_partition_size, silent);

			// Começa do final e volta procurando
			// a lista com mais elementos a serem
			// combinados
			int next = number_of_lists - 1;

			while (next >= 0 &&
					(static_cast<std::vector<int>::size_type>(indices[next] + 1) >= lists_of_attribs[next].size()))
				next--;

			// Nenhuma lista encontrada
			// então não há mais combinações
			if (next < 0)
				break;

			// Se encontrou move-se para o próximo
			// elemento na lista
			indices[next]++;

			// Para todas as listas à direita desta
			// o índice de combinações novamnete retorna
			// para o primeiro elemento
			for (int i = next + 1; i < number_of_lists; i++)
				indices[i] = 0;

		}
	}

	}
}


/**
 * Computa um cubo de dados - essa é uma função fixa do framework
 *
 * @param cube_table o nome de arquivo para o dataset
 * @param num_dims número de dimensões do dataset
 * @param num_meas número de medidas do dataset
 * @param tuple_partition_size o tamanho da partição (por tuplas)
 * @param dim_partition_size o tamanho da partição (por dimensões)
 * @param reading_rate a taxa de leitura do disco usada pelos processos
 * @param tbloc variável específica do algoritmo bCubing
 * @param output_folder arquivo de saída, se o cubo for escrito em disco
 * @param my_rank o rank MPI do processo atual
 * @param queries as consultas que serão executadas
 * @param on_demand booleano que determina se a computação será feita "on-demand"eiros.
 *
 */
void BlockCube::ComputeCube(std::string cube_table, int num_dims,
                int num_meas, int tuple_partition_size, int dim_partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand, std::vector<int> tuple_partition_listings, std::vector<int> dim_partition_listings)
{

		//Um handler de arquivos MPI para o dataset completo
		//Esse arquivo DEVE ser composto apenas de números, onde em cada linha
		//o primeiro número é um inteiro que representa um TID,
		//os n números seguintes são inteiros que representam valores de atributo,
		//os n números seguintes são ponto flutuante e representam medidas
		//Separados por vírgula.
        MPI_File data;

        //Uma forma de saber o endereço de um TID numa tupla do dataset
        MPI_Offset tid_offset;

        //Uma forma de saber o endereço das dimensões numa tupla do dataset
        MPI_Offset dims_offset;

        //Uma forma de saber o endereço das medidas numa tupla do dataset
        MPI_Offset meas_offset;

        //Abra o arquivo fornecido pelo usuário e guarde o endereço dele no handler
        MPI_File_open(MPI_COMM_WORLD, cube_table.c_str(), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &data);

        //Conta quantas tuplas foram lidas para o processo de computar o cubo
        //Esse valor começa em zero e deve terminar quando ler todas as tuplas da partição
        int tuples_read = 0;

        //Ler todas as tuplas da partição de uma só vez seria bem ruim em uso de memória principal
        //Criamos um incremento com base na taxa de ingestão desejada, para ler por partes
        int read_increment = tuple_partition_size / reading_rate;

        //Ajusta o tamanho do buffer de leitura para o tamanho do incremento
        read_buffer.resize(read_increment);

        //Ajusta o tamanho das tuplas no buffer de leitura, que a princípio não se sabe
        //Assim é possível ler tuplas com a quantidade correta de dimensões e medidas
        for (TupleType &tuple : read_buffer)
        {
        	tuple.dims.resize(num_dims);
        	tuple.meas.resize(num_meas);
        }

        ///////////////////// ESPECÍFICO DO ALGORITMO BCUBING /////////////////////

        //ID do primeiro bloco, blocos subsequentes terão IDs incrementais
        int bid = 0;

        //Quantidade de tuplas no bloco identificado por BID, deve ser igual a tbloc
        int bid_tuples_count = 0;

        //Sabendo o tamanho da partição e o tamanho de bloco conseguiremos saber quantos blocos serão criados
        int num_bids = (tuple_partition_size % tbloc == 0) ? tuple_partition_size / tbloc : tuple_partition_size / tbloc + 1;

        //Salva o tamanho do bloco numa variável privada do objeto da classe cubo
        this->tbloc = tbloc;

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBloc
        bbloc.resize(num_dims);

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBlocRAM
        bblocRAM.resize(num_dims);

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingCount (cache de COUNT em memória)
        //Essa estrutura é parte da bCubingBlocRAM mas na implementação achei mais fácil fazer separada
        bblocRAMFreq.resize(num_dims);

        //Garantiamos a alocação de memória para todos os BIDS nessa estrutura auxiliar
        for(auto & dim_bids : bblocRAMFreq){
        	dim_bids.resize(num_bids);
        }

        //Agora que sabemos a quantidade de blocos, redimensiona a estrutura bCubingMeasures
        bmeas.resize(num_bids);

        //Calculamos um bloco por vez apenas, em todas as dimensões, assim utiliza menos memória
        //Depois de processar o bloco será salvo em disco imediatamente
        for (auto & bids_per_dimension : bbloc) {
        	bids_per_dimension.resize(1);
        }

    	//Cada processo MPI tem um diretório próprio para armazenar dados do cubo
    	std::string process_directory = output_folder + "/proc" + std::to_string(my_rank);

        //Se o diretório de saída estiver vazio, significa que foi recém criado (cubo ainda não computado)
        //Nesse caso procede à computação do cubo normalmente
        if(!boost::filesystem::exists(process_directory)){

            //Cria o diretório associado ao processo
            boost::filesystem::create_directory(process_directory);

            //No diretório do processo cria diretórios (inicialmente vazios) para cada uma das dimensões
            //Os blocos usados pelo algoritmo bcubing serão salvos nesses diretórios
            for (int dim_number = 0; dim_number < num_dims; dim_number++)
            {
                    boost::filesystem::create_directory(process_directory + "/dim" + std::to_string(dim_number));
            }

            //Enquanto não ler todas as tuplas dessa partição
            while (tuples_read < tuple_partition_size)
            {
            		//Verifica se na próxima leitura que vai fazer não irá passar do tamanho máximo da partição
            		//Se necessário, reduz o incremento e o buffer de leitura
                    if ((tuples_read + read_increment) > tuple_partition_size)
                    {
                    		//Recalcula o tamanho do incremento apenas com as tuplas restantes
                            read_increment = tuple_partition_size - tuples_read;

                            //Redimensiona, reduzindo potencialmente o tamanho do buffer de leitura
                            read_buffer.resize(read_increment);
                    }

                    //Lê as tuplas do incremento atual e salva no buffer de leitura
                    for (int next_tuple = 0; next_tuple < read_increment; ++next_tuple)
                    {
                    		//Soma o tamanho das partições de todos os os processos anteriores ao atual
                    		//Assim consegue saber em que ponto do arquivo deve começar a ler os dados
                    		int sum_of_partitions = std::accumulate(tuple_partition_listings.begin(), tuple_partition_listings.begin() + my_rank, 0);

                    		//A posição do próximo TID no arquivo binário, indica o começo da próxima tupla que deve ser lida
                    		//Deve pular as tuplas de partições de outros processos e também tuplas que já foram lidas
                    		//As tuplas já lidas devem ser contadas tanto as lidas nesse incremento quanto o total de lidas
                    		//Calcula com base no tamanho em bytes de uma tupla
                            tid_offset =
                                            (sum_of_partitions + (next_tuple + tuples_read))
                                                            * (sizeof(int)
                                                                            + (num_dims
                                                                                            * sizeof(int))
                                                                            + (num_meas
                                                                                            * sizeof(float)));
                            //A posição das próximas dimensões no arquivo binário
                            //Se sabemos a posição do TID, e ele é inteiro, as dimensões começam nessa posição acrescida to tamanho do TID
                            dims_offset = tid_offset + sizeof(int);

                            //A posição das próximas medidas no arquivo binário
                            //Se sabemos a posição das dimensões, e todas são inteiras, as medidas começam nessa posição acrescida dos tamanhos das dimensões
                            meas_offset = dims_offset
                                            + ((num_dims) * sizeof(int));

                            //Leia um inteiro para o TID e guarde no buffer
                            MPI_File_read_at(data, tid_offset, &read_buffer[next_tuple], 1,
                            MPI_INT,
                            MPI_STATUS_IGNORE);

                            //Leia todas as dimensões da tupla e guarde no buffer
                            MPI_File_read_at(data, dims_offset,
                                            &read_buffer[next_tuple].dims[0], 1,
                                            MPI_TUPLE_DIMS, MPI_STATUS_IGNORE);

                            //Leia todas as medidas da tupla e guarde no buffer
                            MPI_File_read_at(data, meas_offset,
                                            &read_buffer[next_tuple].meas[0], 1,
                                            MPI_TUPLE_MEAS,
                                            MPI_STATUS_IGNORE);

                    }

                    //Nesse ponto o buffer de leitura foi preenchido
                    //NÃO leu todas as tuplas da partição, apenas um incremento
                    //Ex: pode ser que o tamanho da partição seja 1 milhão e aqui leu as primeiras 250.000

                    //Para cada uma das tuplas do buffer
                    for (TupleType &tuple : read_buffer)
                    {
                    		//A princípio considere que a tupla é "útil"
                            bool use_tuple = true;

                            //Se cubo sob demanda, verifique se a tupla é realmente "útil"
                            if (on_demand)
                            {
                            		//Quantidade de consultas que envolvem a tupla
                                    int useful_in_queries = 0;

                                    //Para cada consulta a ser executada
                                    for (std::vector<int> &query : queries)
                                    {
                                    		//Se for útil para pelo menos uma consulta é suficiente
                                    		//Não precisa verificar se é util para todas
                                            if (useful_in_queries == 0)
                                            {
                                            		//Valores úteis na tupla para a consulta
                                                    int useful_values = 0;

                                                    //Veja se cada valor da tupla bate com a consulta
                                                    //Também valida se consulta tem agregação ou inquire
                                                    for (int dim_number = 0; dim_number < num_dims;
                                                    		dim_number++)
                                                    {
                                                    		//A tupla é útil se o o valor na dimensão for
                                                    		//exatamente o valor solicitado na consulta ou
                                                    		//se para a dimensão a consulta envolve agregação
                                                    		//ou inquire
                                                            if ((tuple.dims[dim_number] == query[dim_number])
                                                                            || query[dim_number]
                                                                                            == -1
                                                                            || query[dim_number]
                                                                                            == -2)
                                                            {
                                                            		//Encontrou um valor útil na consulta
                                                                    useful_values++;
                                                            }
                                                    }

                                                    //Se todos os valores atenderem a consulta, então a tupla é útil
                                                    if (useful_values == num_dims)
                                                    {
                                                            useful_in_queries++;
                                                    }
                                            }
                                    }

                                    //Se não for útil para nenhuma consulta, marque para não usar
                                    if(useful_in_queries == 0){
                                    	use_tuple = false;
                                    }
                            }

                            //Tupla será inserida no cubo de dados
                            if (use_tuple)
                            {
                            	//Conta uma nova tupla no bloco atual
                            	bid_tuples_count++;

                            	//Se essa nova tupla significar que irá exceder tamanho do bloco
                            	//Crie um novo BID e coloque a tupla nele
                            	if(bid_tuples_count > tbloc){

                            		//Primeiro salva os dados do BID já processado, que estava em memória, para disco
                                    for (int dim_number = 0; dim_number < num_dims; dim_number++)
                                    {
                                    		//SALVA O BLOCO EM DISCO

											//Esse é o diretório onde iremos salvar os dados do bloco recém criado, tanto TIDs quanto medidas
                                    		std::string bloc_directory = process_directory + "/dim" + std::to_string(dim_number) + "/bloc" + std::to_string(bid);

                                    		//Como o bloco acabou de ser computado, não existe em disco ainda, então criamos o diretório
            								boost::filesystem::create_directory(bloc_directory);

            								//Para organizar melhor temos uma pasta específico para os dados de TIDs do bloco
            								boost::filesystem::create_directory(bloc_directory + "/tids/");

            								//Temos também uma pasta específica para os dados de medidas do bloco
            								boost::filesystem::create_directory(bloc_directory + "/meas/");


                                    		//Nome do arquivo onde o BID será salvo - diretório do processo, num diretório específico da dimensão
                                            std::string bid_filename = bloc_directory + "/tids/" + std::to_string(bid) + ".bin";

                                            //Indica que o arquivo de saída será um stream de dados binários
                                            std::ofstream ofb(bid_filename.c_str(), std::ofstream::binary);

                                            //Serialização é feita pela biblioteca Boost
                                            boost::archive::binary_oarchive ob(ofb, boost::archive::no_header);

                                            //Salva os dados do BID recém criado em disco
                                            ob & bbloc[dim_number][0];

                                            //SALVA AS MEDIDAS DO BLOCO EM DISCO

                                            //Nome do arquivo onde as medidas do BID serão salvas - diretório do processo, num diretório específico da dimensão
                                            std::string bid_meas_filename = bloc_directory + "/meas/" + std::to_string(bid) + ".bin";

                                            //Indica que o arquivo de saída será um stream de dados binários
                                            std::ofstream ofm(bid_meas_filename.c_str(), std::ofstream::binary);

                                            //Serialização é feita pela biblioteca Boost
                                            boost::archive::binary_oarchive om(ofm, boost::archive::no_header);

                                            //Salva os dados de medidas do BID recém criado em disco
                                            om & bmeas[0];

                                    }

                                    //Limpa a memória da variável bbloc
                                    bbloc.clear();

                                    //Redimensiona novamente com base na quantidade de dimensões
                                    bbloc.resize(num_dims);

                                    //Limpa a memória da variável bmeas
                                    bmeas.clear();

                                    //Redimensiona novamente com base na quantidade de blocos
                                    bmeas.resize(num_bids);

                                    //Calculamos um bloco por vez apenas, em todas as dimensões, assim utiliza menos memória
                                    //Depois de processar o bloco será salvo em disco imediatamente
                                    for (auto & bids_per_dimension : bbloc) {
                                    	bids_per_dimension.resize(1);
                                    }

                                    //Identificador do novo bloco
                            		bid++;

                            		//Novo bloco já começa com uma tupla, que não foi possível inserir no bloco anterior
                            		bid_tuples_count = 1;
                            	}


                            	//Atualiza as entradas das variáveis bbloc e bblocRAM
                            	//Inclui nessas estruturas os dados da tupla
    							for (int dim_number = 0; dim_number < num_dims; ++dim_number)
    							{
    									//Em cada dimensão, no único bloco existente na posição zero (processamos um bloco por vez),
    									//insere ou atualiza uma entrada com o valor de atributo da dimensão e um novo TID associado
    									bbloc[dim_number][0][tuple.dims[dim_number]].push_back(tuple.tid);

    									//Em cada dimensão insere ou atualiza uma entrada com o valor de atributo da dimensão e um novo BID associado
    									bblocRAM[dim_number][tuple.dims[dim_number]].insert(bid);

										//Atualiza a contagem de frequência do valor de atributo no BID, que fica em memória
    									bblocRAMFreq[dim_number][bid][tuple.dims[dim_number]]++;
    							}

    							//Atualiza as entradas da variável bCubingMeasures
    							//Inclui nessa estrutura os dados da tupla
    							for (int meas_number = 0; meas_number < num_meas; ++meas_number)
    							{
    									//Insere ou atualiza uma entrada com o valor do TID e as medidas associadas
    									//Preenche um bloco por vez na memória e depois salva em disco, por isso aqui sempre usa a primeira posição
    									bmeas[0][tuple.tid].push_back(tuple.meas[meas_number]);
    							}

                            }

                    }

                    //A última escrita de dados da memória pro disco
                    //Ao longo do laço anterior deve ter escrito várias vezes
                    //Aqui é o *último bloco*
                    for (int dim_number = 0; dim_number < num_dims; dim_number++)
                    {
                		//SALVA O BLOCO EM DISCO

						//Esse é o diretório onde iremos salvar os dados do bloco recém criado, tanto TIDs quanto medidas
                		std::string bloc_directory = process_directory + "/dim" + std::to_string(dim_number) + "/bloc" + std::to_string(bid);

                		//Como o bloco acabou de ser computado, não existe em disco ainda, então criamos o diretório
						boost::filesystem::create_directory(bloc_directory);

						//Para organizar melhor temos uma pasta específico para os dados de TIDs do bloco
						boost::filesystem::create_directory(bloc_directory + "/tids/");

						//Temos também uma pasta específica para os dados de medidas do bloco
						boost::filesystem::create_directory(bloc_directory + "/meas/");

                		//Nome do arquivo onde o BID será salvo - diretório do processo, num diretório específico da dimensão
                        std::string bid_filename = bloc_directory + "/tids/" + std::to_string(bid) + ".bin";

                        //Indica que o arquivo de saída será um stream de dados binários
                        std::ofstream ofb(bid_filename.c_str(), std::ofstream::binary);

                        //Serialização é feita pela biblioteca Boost
                        boost::archive::binary_oarchive ob(ofb, boost::archive::no_header);

                        //Salva os dados do BID recém criado em disco
                        ob & bbloc[dim_number][0];

                        //SALVA AS MEDIDAS DO BLOCO EM DISCO

                        //Nome do arquivo onde as medidas do BID serão salvas - diretório do processo, num diretório específico da dimensão
                        std::string bid_meas_filename = bloc_directory + "/meas/" + std::to_string(bid) + ".bin";

                        //Indica que o arquivo de saída será um stream de dados binários
                        std::ofstream ofm(bid_meas_filename.c_str(), std::ofstream::binary);

                        //Serialização é feita pela biblioteca Boost
                        boost::archive::binary_oarchive om(ofm, boost::archive::no_header);

                        //Salva os dados de medidas do BID recém criado em disco
                        om & bmeas[0];

                    }

                    //Pode ser que o incremento de leitura tenha sido alterado acima, então volta para o valor padrão
                    read_increment = tuple_partition_size / reading_rate;

                    //Adiciona a quantidade de tuplas recém lidas ao total de tuplas lidas
                    tuples_read += read_increment;
            }

    		//SALVA INDICE DE BLOCOS DA RAM EM DISCO

            //Nome do arquivo onde o índice atualmente em memória principal será salvo
            std::string bbloc_ram_filename = process_directory + "/" + "index.ram";

            //Indica que o arquivo de saída será um stream de dados binários
            std::ofstream ofram(bbloc_ram_filename.c_str(), std::ofstream::binary);

            //Serialização é feita pela biblioteca Boost
            boost::archive::binary_oarchive oram(ofram, boost::archive::no_header);

            //Salva os dados de medidas do BID recém criado em disco
            oram & bblocRAM;

    		//SALVA CACHE DE FREQUENCIA DA RAM EM DISCO

            //Nome do arquivo onde o índice atualmente em memória principal será salvo
            std::string bbloc_ram_freq_filename = process_directory + "/" + "freq.ram";

            //Indica que o arquivo de saída será um stream de dados binários
            std::ofstream offreq(bbloc_ram_freq_filename.c_str(), std::ofstream::binary);

            //Serialização é feita pela biblioteca Boost
            boost::archive::binary_oarchive ofreq(offreq, boost::archive::no_header);

            //Salva os dados de medidas do BID recém criado em disco
            ofreq & bblocRAMFreq;

        } else { //Cubo já computado

            //Nome do arquivo onde o índice atualmente em memória principal está salvo
			std::string bbloc_ram_filename = process_directory + "/" + "index.ram";

            //Indica que o arquivo de entrada será um stream de dados binários
			std::ifstream ifram(bbloc_ram_filename.c_str(), std::ifstream::binary);

            //Serialização é feita pela biblioteca Boost
			boost::archive::binary_iarchive iram(ifram, boost::archive::no_header);

			//Lê os dados do disco e carrega na variável em memória
			iram & bblocRAM;

            //Nome do arquivo onde o índice atualmente em memória principal está salvo
			std::string bbloc_ram_freq_filename = process_directory + "/" + "freq.ram";

            //Indica que o arquivo de entrada será um stream de dados binários
			std::ifstream iffreq(bbloc_ram_freq_filename.c_str(), std::ifstream::binary);

            //Serialização é feita pela biblioteca Boost
			boost::archive::binary_iarchive ifreq(iffreq, boost::archive::no_header);

			//Lê os dados do disco e carrega na variável em memória
			ifreq & bblocRAMFreq;

        }

        //Limpa o espaço em memória da bCubingBloc
        bbloc.clear();

        //Redimensiona pela quantidade de dimensões pois será usado nas consultas
        bbloc.resize(num_dims);

        //Limpa a memória da variável bmeas
        bmeas.clear();

        //Redimensiona novamente com base na quantidade de blocos
        bmeas.resize(num_bids);

        //Limpa o espaço em memória para o buffer de leitura
        read_buffer.clear();

        //Barreira de sincronismo - só continua quando todos os processos chegarem até aqui
        MPI_Barrier(MPI_COMM_WORLD);

        //Todos os processos fecham a conexão com o arquivo de entrada
        MPI_File_close(&data);
}
