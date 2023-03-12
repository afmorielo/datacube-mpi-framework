#include "DataCube.h"

#ifndef FragCube_h
#define FragCube_h

#include <unordered_map>
#include <set>

class FragCube : public DataCube {
public:
        FragCube(void);
        FragCube(const FragCube& from);
        ~FragCube(void);
        FragCube& operator=(const FragCube& from);
        virtual void ComputeCube(std::string cube_table, int num_dims,
                int num_meas, int tuple_partition_size, int dim_partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand, std::vector<int> tuple_partition_listings, std::vector<int> dim_partition_listings);
        virtual void QueryCube(std::vector<int> query, std::string queries_ops, int my_rank, int num_dims, std::string output_folder, int num_procs);
protected:
private:
        std::vector<std::unordered_map<int, std::vector<int>>> iindex;
        std::unordered_map<int, std::vector<float>> imeas;
        std::map<std::vector<int>,std::map<std::vector<int>, std::vector<int>>> cuboids;
};

#endif // FragCube_h
