#include "FragCube.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <chrono>

void FragCube::PointQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){

}

void FragCube::InquireQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){

}

void FragCube::InquirePointQuery(std::vector<std::vector<int> >& arr, int my_rank, int num_dims, std::string output_folder, int num_procs){

}

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

void FragCube::QueryCube(std::vector<int> query, int my_rank, int num_dims, std::string output_folder, int num_procs){

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
			points_attrs.push_back(query[operand]);
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
		std::vector<int> tids_intersection = cuboids[points][points_attrs];

		//Atualiza o valor de COUNT localmente para a consulta
		local_count = tids_intersection.size();

	} else {
		//Para cada uma das dimensões associadas à consultas point
		for(auto & dim_number : points){
			//Busca os TIDs da estrutura de índice invertido mantida em memória
			lists_of_tids.push_back(cuboids[{ dim_number }][{ points_attrs[dim_number] }]);
		}

		//Faz a interseção das listas de TIDS encontradas
		std::vector<int> tids_intersection = IntersectMultipleVectors(lists_of_tids);

		//Atualiza o valor de COUNT localmente para a consulta
		local_count = tids_intersection.size();
	}

	//Se não tiver nenhuma consulta inquire, deve ter apenas points
	if(inquires.empty()){
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
        std::vector<T>::size_type fragment_size = 2;

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
