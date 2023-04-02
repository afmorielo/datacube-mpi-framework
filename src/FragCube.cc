#include "FragCube.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sstream>

FragCube::FragCube()
{
}

FragCube::FragCube(const FragCube&)
{
}

FragCube::~FragCube()
{
}

FragCube&
FragCube::operator=(const FragCube&)
{
        return *this;
}

void FragCube::QueryCube(std::vector<int> query, std::string queries_ops, std::map<std::vector<int>, int>& query_cache, int my_rank, int num_dims, int num_tuples, std::string output_folder, int num_procs, int tuple_partition_size){

	//Verifica no cache de consultas se a consulta já foi repetida para evitar retrabalho - útil nas inquires
	if(query_cache.count(query) == 0){

	//Listas de TIDs fazem parte da resposta da consulta
	//Cada posição do vector armazena uma lista de TIDs
	std::vector<std::vector<int>> lists_of_tids;

	//Lista final com a interseção de todos os TIDs que respondem à consulta
	std::vector<int> tids_intersection;

	//No caso do inquire irá armazenar listas de valores de atributo
	//Essas listas serão usadas para gerar as consultas pontuais derivadas do inquire
	std::vector<std::vector<int>> lists_of_attribs(num_dims);

	//Lista de dimensões sobre as quais incide um operador inquire
	std::vector<int> inquires;

	//Lista de dimensões sobre as quais incide um operador point
	std::vector<int> points;

	//Lista de atributos buscados sobre os quais incide um operador point
	std::vector<int> points_attrs;

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
			//Atualiza a lista de dimensões afetada por points
			points.push_back(dim);

			//Atualiza a lista de atributos com o valor dessa point (será usado para buscar nos cuboides)
			points_attrs.push_back(query[dim]);
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

	//Normalmente o ideal seria fazer uma busca melhor aqui
	//Verificamos se existe um cuboide que engloba TODAS as dimensões associadas a points
	//No entanto, se não houver o cuboide de mais alto nível poderia combinar cubois menores
	//E.g: a consulta tem points nas dimensões [0,1,3] e tentamos buscar o cuboide [0,1,3], que não existe
	//o ideal seria então buscar os cuboides [0,1],[0,3] e [1,3] e só em último caso os cuboides [0],[1],[3]
	//Mas aqui tentamos o cuboide de mais alto nível e depois partimos direto pro mais baixo nível
	if(cuboids.count(points) != 0){
		//Busca a lista de TIDs que já foi pré-calculada
		tids_intersection = cuboids[points][points_attrs];

		//Atualiza o valor de COUNT localmente para a consulta
		local_count = tids_intersection.size();

	} else {

		//Para cada uma das dimensões associadas à consultas point buscamos nos cuboides
		//a lista de TIDs associada ao atributo requerido na consulta para essa dimensão
		for(std::vector<T>::size_type index = 0; index != points.size(); index++) {
			//Busca os TIDs da estrutura de índice invertido mantida em memória
			lists_of_tids.push_back(cuboids[{ points[index] }][{ points_attrs[index] }]);
		}

		//Faz a interseção das listas de TIDS encontradas
		tids_intersection = IntersectMultipleVectors(lists_of_tids);

		//Atualiza o valor de COUNT localmente para a consulta
		local_count = tids_intersection.size();
	}

	//Se não tiver nenhuma consulta inquire, deve ter apenas points
	if(inquires.empty()){

		//Somente será executado quando todos processos chegarem nesse ponto
		//Combina as respostas de todos, e todos sabem qual a contagem global de tuplas
		MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//Se essa consulta teve algum resultado, vamos guardá-la em cache para evitar a reconsulta
		if(global_count > 0){
			query_cache[query] = global_count;
		}

		//Apresenta inicialmente a contagem de tuplas
		if(my_rank == 0){
			if(global_count > 0){
				for(auto & operand : query) {
					(operand == -1) ? std::cout << '*' << ' ' : std::cout << operand << ' ';
				}
				std::cout << ": " << global_count << ' ';
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

							//Acessa o índice em memória para obter o valor de medida e faz a SOMA
							local_sum += imeas[tid][measure_id];

						}

	        		}

	        		//Somente será executado quando todos processos chegarem nesse ponto
	        		//Combina as respostas de todos num só processo
	        		MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	        		//Apresenta o valor agregado de SOMA das medidas
	        		if(my_rank == 0){
	        			if(global_count > 0){
	        				std::cout << "Sm(M" << measure_id << ") = " << global_sum << ' ';
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

							//Obtém a lista de medidas localmente
							local_measures.push_back(imeas[tid][measure_id]);

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

	        			    //Obtemos o valor mediano, após a ordenação
	        				std::cout << "Md(M" << measure_id << ") = " << global_measures[n] << ' ';
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
				//Oficialmente terminamos de responder a consulta!
				std::cout << std::endl;
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

		//Para cada uma das dimensões associadas à consultas inquire
		for(auto & dim_number : inquires){

			//Irá guardar os valores de atributo possíveis na dimensão
			std::set<int> attribs_set;

			//Um valor extra necessário é o de agregação '*' pois não está nos dados e é usado em consulta válida de inquire
			//E.g: uma dimensão com atributos [1,2] na verdade tem os atributos [-1,1,2], onde -1 indica TODOS
			attribs_set.insert(-1);

			//Para cada par associando um valor de atributo a uma lista de TIDs
			for (auto& attrib_pair_tids: cuboids[{ dim_number }]) {
				//Extraia apenas o valor de atributo, o primeiro elemento do pair, e salve na lista de atributos
				attribs_set.insert(attrib_pair_tids.first[0]);
			}

			//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
			std::vector<int> attribs_vector(attribs_set.begin(), attribs_set.end());

			//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
			lists_of_attribs[dim_number] = attribs_vector;

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
			QueryCube(query, queries_ops, query_cache, my_rank, num_dims, num_tuples, output_folder, num_procs, tuple_partition_size);

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

void FragCube::ComputeCube(std::string cube_table, int num_dims,
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

        ///////////////////// ESPECÍFICO DO ALGORITMO FRAGCUBING /////////////////////

        //Tamanho de fragmento, ou seja, quantas dimensões serão usadas em cada cuboide
        std::vector<T>::size_type fragment_size = 1;

        //Dimensões de um determinado cuboide, se a divisão for exata deve ter sempre o tamanho do fragmento definido
        std::vector<int> fragment_dims;

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura Inverted Index
        iindex.resize(num_dims);

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


                        if (use_tuple)
                        {
                                for (int dim_number = 0; dim_number < num_dims; ++dim_number)
                                {
                                        iindex[dim_number][tuple.dims[dim_number]].push_back(tuple.tid);
                                }

                                for (int meas_number = 0; meas_number < num_meas; ++meas_number)
                                {
                                        imeas[tuple.tid].push_back(tuple.meas[meas_number]);
                                }
                        }

                }

                //Pode ser que o incremento de leitura tenha sido alterado acima, então volta para o valor padrão
                read_increment = tuple_partition_size / reading_rate;

                //Adiciona a quantidade de tuplas recém lidas ao total de tuplas lidas
                tuples_read += read_increment;

        }

        //Agora, a partir dos índices, pode pré-computar cuboides
        //Passamos novamente por todas as dimensões para isolar aquelas que fazem parte de cada cuboide
        for (int dim_number = 0; dim_number < num_dims; ++dim_number){

        	//A dimensão é incluída em um fragmento
        	fragment_dims.push_back(dim_number);

        	//Verificamos - o fragmento atingiu o tamanho desejado (ou acabaram as dimensões e deverá ser um fragmento um pouco menor?)
        	if(fragment_dims.size() == fragment_size || dim_number == num_dims - 1){

				//Aqui isolamos as dimensões que fazem parte desse fragmento
				//A quantidade de cuboides para esse fragmento é exponencial com base no seu tamanho
				int fragment_cuboids_total = IntegerPow(2, fragment_dims.size());

				//Vamos calcular todos os cuboides desse fragmento
				for(int cuboid_number = 1 ; cuboid_number < fragment_cuboids_total ; cuboid_number++){

					//Aqui obtemos o numero de cuboide, mas na verdade isso será tratado como inteiro (e.g: 00000001)
					int binary = cuboid_number;

					//Essa é a lista de dimensões desse cuboide, isso varia, por exemplo (A,B,AB,AC,BC,ABC)
					std::vector<int> cuboid_dims;

					//Aqui fazemos uma contagem binária para gerar as combinações com um mínimo de memória e tempo
					//00000001 = A
					//00000010 = B
					//00000011 = AB
					//etc.
					for(std::vector<T>::size_type dim_index = 0 ; dim_index < fragment_dims.size() ; dim_index++){

						//Se o bit for 1 então será usado na combinação
						if (binary & 1){    // Verifica um bit do número.
							cuboid_dims.push_back(fragment_dims[dim_index]);
						}

						binary = binary >> 1;  //Faz a rotação de bits para ir para a próxima combinação
					}

					//Gerar o cuboide é equivalente a executar inquires nas dimensões, então obtemos valores de atributo de cada dimensão
					//Essas listas serão usadas para gerar as consultas pontuais derivadas do inquire
					std::vector<std::vector<int>> lists_of_attribs;

					//Para cada uma das dimensões associadas ao cuboide
					for(auto & dim_number : cuboid_dims){

						//Irá guardar os valores de atributo possíveis na dimensão
						std::set<int> attribs_set;

						//NOS CUBOIDES NÃO É INCLUIDO O VALOR DE AGREGAÇÃO
						//Um valor extra necessário é o de agregação '*' pois não está nos dados e é usado em consulta válida de inquire
						//E.g: uma dimensão com atributos [1,2] na verdade tem os atributos [-1,1,2], onde -1 indica TODOS
						//attribs_set.insert(-1);

						//Para cada par associando um valor de atributo a uma lista de TIDs
						for (auto& attrib_pair_tids: iindex[dim_number]) {
							//Extraia apenas o valor de atributo, o primeiro elemento do pair, e salve na lista de atributos
							attribs_set.insert(attrib_pair_tids.first);
						}

						//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
						std::vector<int> attribs_vector(attribs_set.begin(), attribs_set.end());

						//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
						lists_of_attribs.push_back(attribs_vector);

					}

					//Agora serão geradas as combinações de atributos para o cuboide, equivalente ao inquire
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

					while (1) {

						//Listas de TIDs fazem parte da celula do cuboide, então precisaremos determiná-las
						//Cada posição do vector armazena uma lista de TIDs
						std::vector<std::vector<int>> lists_of_tids;

						//Armazena a celula do cuboide gerada pelo processo nessa iteração
						std::vector<int> cuboid_cell(lists_of_attribs.size());

						//Gera a consulta com base na próxima combinaçao de atributos
						for (int i = 0; i < number_of_lists; i++){
							cuboid_cell[i] = lists_of_attribs[i][indices[i]];
						}

						//Agora que temos os valores de atributo da célula, iremos buscar as listas de TIDs
						for(std::vector<T>::size_type attr_index = 0; attr_index < cuboid_cell.size(); attr_index++){

							//Busca os TIDs da estrutura de índice invertido mantida em memória
							lists_of_tids.push_back(iindex[cuboid_dims[attr_index]][cuboid_cell[attr_index]]);

						}

						//Faz a interseção das listas de TIDS encontradas
						std::vector<int> tids_intersection = IntersectMultipleVectors(lists_of_tids);

						//Somente são inseridos no cuboide quando a interseção não é vazia
						if(tids_intersection.size() > 0){
							cuboids[cuboid_dims][cuboid_cell] = tids_intersection;
						}

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

        		//Limpa o espaço em memória para que seja usado pelo próximo fragmento
        		fragment_dims.clear();

        	}

        }

        //Limpa o espaço em memória do índice invertido pois agora já computou os cuboides
        iindex.clear();

        //Limpa o espaço em memória para o buffer de leitura
        read_buffer.clear();

        //Barreira de sincronismo - só continua quando todos os processos chegarem até aqui
        MPI_Barrier(MPI_COMM_WORLD);

        //Todos os processos fecham a conexão com o arquivo de entrada
        MPI_File_close(&data);

}
