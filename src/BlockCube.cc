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

void BlockCube::QueryCube(std::vector<int> query, std::string queries_ops, int my_rank, int num_dims, std::string output_folder, int num_procs){

	//Cada processo MPI tem um diretório próprio para armazenar dados do cubo
	std::string process_directory = output_folder + "/proc" + std::to_string(my_rank);

	//Listas de BIDs associados a valores de atributo da consulta
	//Cada posição do vector armazena uma lista de BIDs
	std::vector<std::vector<int>> lists_of_bids;

	//Listas de TIDs fazem parte da resposta da consulta
	//Cada posição do vector armazena uma lista de TIDs
	std::vector<std::vector<int>> lists_of_tids;

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

	//Faz a interseção das listas de BIDS encontradas
	std::vector<int> bids_intersection = IntersectMultipleVectors(lists_of_bids);

	//Se não tiver nenhuma consulta inquire, deve ter apenas points
	if(inquires.empty()){

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
				}

				//Vai armazenar uma quantidade de listas igual à quantidade de consultas do tipo point
				lists_of_tids.push_back(tids);

			}

			//Faz a interseção das listas de TIDS encontradas
			std::vector<int> tids_intersection = IntersectMultipleVectors(lists_of_tids);

			//Atualiza o valor de COUNT localmente para a consulta
			local_count = tids_intersection.size();
		}

		//Somente será executado quando todos processos chegarem nesse ponto
		//Combina as respostas de todos num só processo
		MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if(my_rank == 0){
			if(global_count > 0){
				for(auto & operand : query) {
					(operand == -1) ? std::cout << '*' << ' ' : std::cout << operand << ' ';
				}
				std::cout << ": " << global_count << std::endl;
			}
		}
	} else { //Consulta tem pelo menos um operador inquire

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
				}

				//Converte o conjunto para um vector simples, agora que sabemos que está ordenado (para poder ser usado no método de interseção)
				std::vector<int> attribs_vector(attribs_set.begin(), attribs_set.end());

				//Insere a lista de atributos da dimensão com inquire na posição associada à dimensão nas listas de atributos
				lists_of_attribs[dim_number] = attribs_vector;

			}


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

		//Quantidade de processos que já terminaram de gerar consultas
		//Inicialmente é zero pois nenhuma consulta foi gerada/respondida ainda
		int procs_finished = 0;

		//Armazena uma lista dos indicadores de finalização de cada processo
		//É útil para os demais processos saberem quando alguém vai enviar uma nova consulta ou não
		std::vector<int> procs_finished_control(num_procs);

		//Equivalente a um booleano, 0 - FALSE, 1 -TRUE, que indica se o processo atual terminou de gerar consultas
		int finished = 0;

		//Algumas vezes o processo nem gerou as listas pois a interseção de BIDs foi vazia
		//Isso significa que a partição que o processo recebeu não ajuda na consulta
		//Sendo assim, ele já finalizou e verificar a primeira lista é suficiente pra confirmar isso
		if(lists_of_attribs.front().empty()){
			finished = 1;
		}

		//Só sai desse laço com um "break"
		while (1) {

			//Soma os indicadores de finalização de todos os processos e salva numa variável com a contagem geral, todos tem uma cópia
			MPI_Allreduce(&finished, &procs_finished, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			//Se todos tiverem finalizado, todos terão definido um valor "1" e a soma será igual ao número de processos
			//Nesse caso pode encerrar a execução
			if(procs_finished == num_procs)
				break;

			//Captura os indicadores de finalização de todos os processos e salva na lista geral, que todos tem uma cópia
			MPI_Allgather(&finished, 1, MPI_INT, procs_finished_control.data(), 1, MPI_INT, MPI_COMM_WORLD);

			//Armazena a consulta gerada pelo processo nessa iteração
			std::vector<int> query(num_dims);

			//Armazena a consulta, dentre aquelas gerados por todos os processos nessa iteração, sendo executada nesse momento
			std::vector<int> query_running(num_dims);

			//Só faz sentido tentar gerar uma nova consulta se ainda não tiver finalizado
			if(finished == 0){
				//Gera a consulta com base na próxima combinaçao de atributos
				for (int i = 0; i < number_of_lists; i++){
					query[i] = lists_of_attribs[i][indices[i]];
				}
			}

			//Nesse ponto cada processo que poderia gerar uma nova consulta gerou uma nova consulta nessa iteração
			for(int i = 0; i < num_procs; i++){

				//A cada etapa do laço define um processo como o "mandante"
				if(my_rank == i){
					//A consulta que será executada é a consulta gerada por aquele processo
					query_running = query;
				}

				//Todos os processos verificam se o "mandante" já está finalizado
				if(procs_finished_control[i] == 0){

					//A query do "mandante" é enviada para todos os demais processos
					MPI_Bcast(&query_running[0], num_dims, MPI_INT, i, MPI_COMM_WORLD);

					//Todos os processos executam a query
					QueryCube(query_running, queries_ops, my_rank, num_dims, output_folder, num_procs);
				}
			}

			// Começa do final e volta procurando
			// a lista com mais elementos a serem
			// combinados
			int next = number_of_lists - 1;

			if(finished == 0){
				while (next >= 0 &&
					  (static_cast<std::vector<int>::size_type>(indices[next] + 1) >= lists_of_attribs[next].size()))
					next--;

				// Nenhuma lista encontrada
				// então não há mais combinações
				if (next < 0)
					finished = 1;
			}

			if(finished == 0){
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

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBloc
        bbloc.resize(num_dims);

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBlocRAM
        bblocRAM.resize(num_dims);

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
        if(boost::filesystem::is_empty(output_folder)){

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

        } else { //Cubo já computado

            //Nome do arquivo onde o índice atualmente em memória principal está salvo
			std::string bbloc_ram_filename = process_directory + "/" + "index.ram";

            //Indica que o arquivo de entrada será um stream de dados binários
			std::ifstream ifram(bbloc_ram_filename.c_str(), std::ifstream::binary);

            //Serialização é feita pela biblioteca Boost
			boost::archive::binary_iarchive iram(ifram, boost::archive::no_header);

			//Lê os dados do disco e carrega na variável em memória
			iram & bblocRAM;

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
