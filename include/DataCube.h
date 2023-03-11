/*
 * Nesse arquivo temos uma definição de classe base abstrata para um Cubo de Dados
 * A ideia é definir uma classe genérica que tenha os atributos e métodos que qualquer Cubo de Dados precise implementar.
 * Cada algoritmo pode implementar o cubo à sua própria maneira.
 */
#ifndef DataCube_h
#define DataCube_h

#include <vector>
#include <mpi.h>
#include <nadeau.h>

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
        virtual void ComputeCube(std::string cube_table, int num_dims,
                int num_meas, int tuple_partition_size, int dim_partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand, std::vector<int> tuple_partition_listings, std::vector<int> dim_partition_listings) = 0;
        virtual void QueryCube(std::vector<int> query, std::string queries_ops, int my_rank, int num_dims, std::string output_folder, int num_procs) = 0;
        virtual std::vector<int> IntersectTwoVectors(std::vector<int> const& vector_a, std::vector<int> const& vector_b);
        virtual std::vector<int> IntersectMultipleVectors(std::vector<std::vector<int>> &sorted_vectors);
        virtual int IntegerPow(int base, int exp);
        std::vector<TupleType> read_buffer;
        std::vector<int> tuple_partition_listings;
        std::vector<int> dim_partition_listings;
protected:
private:
};

#endif // DataCube_h
