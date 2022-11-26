/*
 * Nesse arquivo temos uma definição de classe base abstrata para um Cubo de Dados
 * A ideia é definir uma classe genérica que tenha os atributos e métodos que qualquer Cubo de Dados precise implementar.
 * Cada algoritmo pode implementar o cubo à sua própria maneira.
 */
#ifndef DataCube_h
#define DataCube_h

#include <vector>
#include <mpi.h>

// Definição de uma tupla como um TID inteiro, uma sequência de valores inteiros para as dimensões
// e uma sequência de números com ponto flutuante para as medidas.
using TupleType =
struct T
{
        int tid;
        std::vector<int> dims;
        std::vector<float> meas;
};

// Esses tipos de dados são permitidos no MPI para facilitar o tratamento
// de estruturas mais complexos do que inteiros e floats simples (variáveis globais)
extern MPI_Datatype MPI_TUPLE_DIMS, MPI_TUPLE_MEAS;

class DataCube
{
public:
        virtual ~DataCube();
        virtual void ComputeCube(std::string cube_table, int num_dims, int num_meas, int partition_size, int ingestion_rate, int tbloc, std::string output_folder, int my_rank,
                        std::vector<std::vector<int>> queries, bool on_demand) = 0;
        virtual void QueryCube(std::vector<std::vector<int>> queries, int my_rank, int num_dims, std::string output_folder, int num_procs) = 0;
        std::vector<TupleType> read_buffer;
protected:
private:
};

#endif // DataCube_h
