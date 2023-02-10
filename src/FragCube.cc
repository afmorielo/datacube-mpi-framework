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
	std::cout << "Consulta" << std::endl;
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
        		if(my_rank == 0){
					std::cout << "Dimensões [";
					for(auto & e : fragment_dims){
						std::cout << e << " ";
					}
					std::cout << "]" << std::endl;

					//A quantidade de cuboides para esse fragmento é exponencial com base no seu tamanho
					int fragment_cuboids = IntegerPow(2, fragment_dims.size());

					//Vamos calcular todos os cuboides desse fragmento
					for(int cuboid_number = 1 ; cuboid_number < fragment_cuboids ; cuboid_number++){

						//Aqui obtemos o numero de cuboide, mas na verdade isso será tratado como inteiro (e.g: 00000001)
						int binary = cuboid_number;

						std::vector<int> cuboid_dims;

						std::cout << "Cuboide >>>> [";

						//Aqui fazemos uma contagem binária para gerar as combinações com um mínimo de memória e tempo
						//00000001 = A
						//00000010 = B
						//00000011 = AB
						//etc.
						for(std::vector<T>::size_type dim_index = 0 ; dim_index < fragment_dims.size() ; dim_index++){

							//Se o bit for 1 então será usado na combinação
							if (binary & 1){    // Verifica um bit do número.
								std::cout << fragment_dims[dim_index] << " ";
								cuboid_dims.push_back(fragment_dims[dim_index]);
							}

							binary = binary >> 1;  //Faz a rotação de bits para ir para a próxima combinação
						}
						std::cout << "]" << std::endl;

						//No caso do inquire irá armazenar listas de valores de atributo
						//Essas listas serão usadas para gerar as consultas pontuais derivadas do inquire
						std::vector<std::vector<int>> lists_of_attribs;

						//Para cada uma das dimensões associadas à consultas inquire
						for(auto & dim_number : cuboid_dims){

							//Irá guardar os valores de atributo possíveis na dimensão
							std::set<int> attribs_set;

							//Um valor extra necessário é o de agregação '*' pois não está nos dados e é usado em consulta válida de inquire
							//E.g: uma dimensão com atributos [1,2] na verdade tem os atributos [-1,1,2], onde -1 indica TODOS
							attribs_set.insert(-1);

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

					    while (1) {

							//Listas de TIDs fazem parte da resposta da consulta
							//Cada posição do vector armazena uma lista de TIDs
							std::vector<std::vector<int>> lists_of_tids;

							//Armazena a consulta gerada pelo processo nessa iteração
							std::vector<int> query(lists_of_attribs.size());

							std::cout << "Celula >>>>>>>> [";

							//Gera a consulta com base na próxima combinaçao de atributos
							for (int i = 0; i < number_of_lists; i++){
								std::cout << lists_of_attribs[i][indices[i]] << " ";
								query[i] = lists_of_attribs[i][indices[i]];
							}

							std::cout << "]: ";

							//Para cada uma das dimensões associadas à consultas point
							for(int attr_index = 0; attr_index < query.size(); attr_index++){

								//Armazena os TIDs associados ao valor de atributo de uma das dimensões da consulta
								std::vector<int> tids;

								std::vector<int> tids_bbloc;

								//Acessa uma única vez a estrutura bCubingBloc para obter os TIDS para esse BID
								if(query[attr_index] != -1){
									tids_bbloc = iindex[cuboid_dims[attr_index]][query[attr_index]];
								} else {
									tids_bbloc.resize(tuple_partition_size);
									std::iota(tids_bbloc.begin(), tids_bbloc.end(), 1);
								}


								//Adiciona às listas de TID os TIDs da dimensão associada à point, para esse BID, com base no valor da consulta
								tids.insert(std::end(tids), std::begin(tids_bbloc), std::end(tids_bbloc));

								//Vai armazenar uma quantidade de listas igual à quantidade de consultas do tipo point
								lists_of_tids.push_back(tids);

							}

							//Faz a interseção das listas de TIDS encontradas
							std::vector<int> tids_intersection = IntersectMultipleVectors(lists_of_tids);

							//Atualiza o valor de COUNT localmente para a consulta
							std::cout << tids_intersection.size() << std::endl;

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

        		fragment_dims.clear();

        	}

        }

        //Limpa o espaço em memória para o buffer de leitura
        read_buffer.clear();

        //Barreira de sincronismo - só continua quando todos os processos chegarem até aqui
        MPI_Barrier(MPI_COMM_WORLD);

        //Todos os processos fecham a conexão com o arquivo de entrada
        MPI_File_close(&data);

}
