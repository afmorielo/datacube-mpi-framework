#ifndef Handler_h
#define Handler_h

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

class Handler {
public:
	bool ParseInput(int argc, char *argv[], int my_rank, int num_procs,
			int &num_dims, int &num_meas, int &num_tuples,
			int &tuple_partition_size, int &dim_partition_size,
			int &reading_rate, int &tbloc, std::string &output_folder,
			std::vector<std::vector<int>> &queries, std::string &cube_algorithm,
			std::string &cube_table, bool &on_demand);
protected:
private:
};

#endif // Handler_h
