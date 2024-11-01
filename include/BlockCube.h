#include "DataCube.h"

#ifndef BlockCube_h
#define BlockCube_h

#include <unordered_map>
#include <set>

class BlockCube: public DataCube
{
public:
        BlockCube(void);
        BlockCube(const BlockCube& from);
        ~BlockCube(void);
        BlockCube& operator=(const BlockCube& from);
        virtual void ComputeCube(std::string cube_table, int num_dims,
                int num_meas, int tuple_partition_size, int dim_partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand, std::vector<int> tuple_partition_listings, std::vector<int> dim_partition_listings);
        virtual void QueryCube(std::vector<int> query, std::string queries_ops, std::map<std::vector<int>, int>& query_cache, int my_rank, int num_dims, int num_tuples, std::string output_folder, int num_procs, int tuple_partition_size, bool silent);
protected:
private:
        std::vector<std::vector<std::unordered_map<int, std::vector<int>>>> bbloc;
        std::vector<std::unordered_map<int, std::set<int>>> bblocRAM;
        std::vector<std::unordered_map<int, std::vector<float>>> bmeas;
        std::vector<std::vector<std::unordered_map<int, int>>> bblocRAMFreq;
        int tbloc = 2; //Tamanho de bloco, o valor padrão é 2
};

#endif // BlockCube_h
