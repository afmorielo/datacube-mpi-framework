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

        //Um objeto da classe Handler, cuja função é garantir que as variáveis passadas pelo usuário sejam válidas
        Handler h;

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
        if (!h.ParseInput(argc, argv, my_rank, num_procs, num_dims, num_meas,
                        num_tuples, tuple_partition_size, dim_partition_size, reading_rate, tbloc, output_folder, queries,
                        cube_algorithm, cube_table, on_demand))
        {
                return 1;
        }

        //Para permitir a leitura das dimensões do disco em segmentos de mesmo tamanho
        MPI_Type_contiguous(num_dims, MPI_INT, &MPI_TUPLE_DIMS);
        MPI_Type_commit(&MPI_TUPLE_DIMS);

        //Para permitir a leitura das tuplas do disco em segmentos de mesmo tamanho
        MPI_Type_contiguous(num_meas, MPI_FLOAT, &MPI_TUPLE_MEAS);
        MPI_Type_commit(&MPI_TUPLE_MEAS);

        //Execução do algoritmo bCubing - computa e depois executa as consultas
        if (cube_algorithm == "bcubing")
        {
                begin_compute = std::chrono::steady_clock::now();
                cube = std::make_unique<BlockCube>(BlockCube());
                cube->ComputeCube(cube_table, num_dims, num_meas, tuple_partition_size, reading_rate, tbloc, output_folder, my_rank, queries, on_demand);
                end_compute = std::chrono::steady_clock::now();

                begin_query = std::chrono::steady_clock::now();
                cube->QueryCube(queries, my_rank, num_dims, output_folder, num_procs);
                end_query = std::chrono::steady_clock::now();
        }
        else if (cube_algorithm == "fragcubing") //Execução do algoritmo fragCubing - computa e depois executa as consultas
        {
                begin_compute = std::chrono::steady_clock::now();
                cube = std::make_unique<FragCube>(FragCube());
                cube->ComputeCube(cube_table, num_dims, num_meas, tuple_partition_size, reading_rate, tbloc, output_folder, my_rank, queries, on_demand);
                end_compute = std::chrono::steady_clock::now();

                begin_query = std::chrono::steady_clock::now();
                cube->QueryCube(queries, my_rank, num_dims, output_folder, num_procs);
                end_query = std::chrono::steady_clock::now();
        }

        if(my_rank==0){
                std::cout << "Time compute = " << std::chrono::duration_cast<std::chrono::seconds> (end_compute - begin_compute).count() << "[s]" << std::endl;
                std::cout << "Time query = " << std::chrono::duration_cast<std::chrono::microseconds> (end_query - begin_query).count() << "[µs]" << std::endl;
        }

        MPI_Type_free(&MPI_TUPLE_DIMS);
        MPI_Type_free(&MPI_TUPLE_MEAS);
        MPI_Finalize();

        return 0;
}
