#include <iostream> 		//Funções de entrada/saída
#include <mpi.h>			//Funções do MPI (utilizo OpenMPI)
#include <chrono>			//Para contagem do tempo de execução
#include "BlockCube.h"		//bCubing - métodos de computação, consulta, atualização
#include "FragCube.h"		//fragCubing - métodos de computação, consulta, atualização
#include "Handler.h"	//Gerenciamento de parâmetros do programa

//Alias para tempo
using Time = std::chrono::steady_clock;

//Alias para um ponto no tempo
using TimePoint = std::chrono::steady_clock::time_point;

//Esse programa precisa de parâmetros por linha de comando!
int main(int argc, char *argv[])
{

		//Um ponteiro inteligente para um cubo de dados (deleta automaticamente ao sair de escopo)
    	std::unique_ptr<DataCube> cube;

		//Armazena o nome do algoritmo que será executado
        std::string cube_algorithm;

        //Armazena o caminho do sistema de arquivos onde a tabela de fatos está
		std::string cube_table;

		//Diretório de saída do cubo, útil para persistir dados em disco
        std::string output_folder;

		//Uma forma de listar todas as consultas a serem executadas no cubo
        std::vector<std::vector<int>> queries;

        //Uma ou mais listas de operações do usuário para executar em medidas cubo de dados
        std::vector<std::string> queries_ops;

        //Esse é um "cache" de consultas respondidas que mapeia a consulta (sequência de inteiros) à quantidade de TIDs encontrados
        //A ideia é que cada consulta que tenha uma frequência de TIDs maior que 0 seja salva aqui para evitar retrabalho
        std::map<std::vector<int>, int> query_cache;

        //Tempos de computação/consultas (para gelar o relatório de execução)
        TimePoint begin_compute; //Começou a computar (ponto no tempo)
        TimePoint end_compute; //Terminou de computar (ponto no tempo)
        TimePoint begin_query; //Começou a consultar (ponto no tempo)
        TimePoint end_query; //Terminou de consultar (ponto no tempo)

        //Se TRUE, o cubo será computado sob demanda. Se FALSE, será computado para qualquer consulta
        bool on_demand;

        int my_rank; 	//Rank do processo (0 até N-1)
        int num_procs;  //Numero de processos em execução
        int num_dims;	//Número de dimensões do dataset/cubo
        int num_meas;	//Número de medidas do dataset/cubo
        int num_tuples;	//Número de tuplas do dataset/cubo
        int tuple_partition_size; //Tamanho da partição, ou seja, divisão das tuplas pelos processos
        int dim_partition_size; //Tamanho da partição, ou seja, divisão das dimensões pelos processos
        int reading_rate;	//Taxa de leitura dos dados de disco para memória
        int tbloc; //Tamanho do bloco, só é aplicável para o algoritmo bCubing

        //Esse são vectors que guardam o tamanho de partição de cada processo (seja por tuplas ou dimensões)
        //Se na posição 0 temos o valor 10, significa que o processo 0 usa uma partiçãod de tamanho 10.
        //Dependendo da entrada fornecida pelo usuário isso signifca 10 tuplas ou 10 dimensões.
        std::vector<int> tuple_partition_listings;
        std::vector<int> dim_partition_listings;

        //Um objeto da classe Handler, cuja função é garantir que as variáveis passadas pelo usuário sejam válidas
        Handler IO;

        //Iniciando ambiente MPI
        MPI_Init(&argc, &argv);

        //Definindo error handler do MPI
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

        //Calculando quantidade de processos em execução
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        //Calculando o rank do processo em execução neste código/instância
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        //Verifica todas as variáveis passadas pelo usuário e garanta que sejam válidas
        //Também atribui valores as variáveis que depois serão usadas para computar e consultar cubos.
        if (!IO.ParseInput(argc, argv, my_rank, num_procs, num_dims, num_meas,
                        num_tuples, tuple_partition_size, dim_partition_size, reading_rate, tbloc, output_folder, queries, queries_ops,
                        cube_algorithm, cube_table, on_demand))
        {
        	MPI_Finalize();
        }
        else{

        	//Aloca a memória necessária para os dados de todos os processos
        	tuple_partition_listings.resize(num_procs);
        	dim_partition_listings.resize(num_procs);

        	//Salva os dados de tamanho de partições, todos os processos tem uma cópia disso
            MPI_Allgather(&tuple_partition_size, 1, MPI_INT, tuple_partition_listings.data(), 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Allgather(&dim_partition_size, 1, MPI_INT, dim_partition_listings.data(), 1, MPI_INT, MPI_COMM_WORLD);

            //Para permitir a leitura das dimensões do disco em segmentos de mesmo tamanho
            MPI_Type_contiguous(num_dims, MPI_INT, &MPI_TUPLE_DIMS);
            MPI_Type_commit(&MPI_TUPLE_DIMS);

            //Para permitir a leitura das tuplas do disco em segmentos de mesmo tamanho
            MPI_Type_contiguous(num_meas, MPI_FLOAT, &MPI_TUPLE_MEAS);
            MPI_Type_commit(&MPI_TUPLE_MEAS);

            //Execução do algoritmo bCubing - computa e depois executa as consultas
            if (cube_algorithm == "bcubing")
            {
                    //Define o ponteiro único para um cubo do tipo desejado
                    cube = std::make_unique<BlockCube>(BlockCube());
            }
            else if (cube_algorithm == "fragcubing") //Execução do algoritmo fragCubing - computa e depois executa as consultas
            {
                	//Define o ponteiro único para um cubo do tipo desejado
                	cube = std::make_unique<FragCube>(FragCube());
            }

            //Essa lista fica salva como parte integrante dos dados de qualquer cubo
            cube->tuple_partition_listings = tuple_partition_listings;

            //Essa lista fica salva como parte integrante dos dados de qualquer cubo
            cube->dim_partition_listings = dim_partition_listings;

    		//Tempo de início da computação
            begin_compute = Time::now();

            //Invoca a implementação do método de computação do cubo
            cube->ComputeCube(cube_table, num_dims, num_meas, tuple_partition_size, dim_partition_size, reading_rate, tbloc, output_folder, my_rank, queries, on_demand, tuple_partition_listings, dim_partition_listings);

    		//Tempo logo após finalização da computação
            end_compute = Time::now();

            //Obtém o uso de memória do processo para o cubo computado
            std::cout << "(" << my_rank << "): " << "Memory usage (bytes) => " << getPeakRSS() << std::endl;

            //O processo principal apresenta dados da computação
            if(my_rank==0){
            		//Espaço deixado intencionalmente em branco
            		std::cout << std::endl;

            		//Apresenta um relatório simples de tempo e uso de memória para computação
                    std::cout << "Cube computed in " << std::chrono::duration_cast<std::chrono::milliseconds> (end_compute - begin_compute).count() << "[ms]" << std::endl;

            		//Espaço deixado intencionalmente em branco
                    std::cout << std::endl;
            }

        	//Irá avaliar consulta a consulta que o usuário forneceu
            for (std::vector<T>::size_type num_query = 0; num_query != queries.size(); num_query++)
            {
        		//Tempo de início das consulta
                begin_query = Time::now();

                //Invoca a implementação do método de consulta do cubo
            	cube->QueryCube(queries[num_query], queries_ops[num_query], query_cache, my_rank, num_dims, num_tuples, output_folder, num_procs, tuple_partition_size);

        		//Tempo logo após finalização da consulta
                end_query = Time::now();

                if(my_rank==0){
                		//Apresenta um relatório simples de tempo de resposta da consulta
                        std::cout << "Query solved in " << std::chrono::duration_cast<std::chrono::microseconds> (end_query - begin_query).count() << "[µs]" << std::endl;

                		//Espaço deixado intencionalmente em branco
                        std::cout << std::endl;
                }
            }

            MPI_Type_free(&MPI_TUPLE_DIMS);
            MPI_Type_free(&MPI_TUPLE_MEAS);
            MPI_Finalize();
        }
}
