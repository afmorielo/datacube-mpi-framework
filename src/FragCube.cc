#include "FragCube.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <chrono>

auto tbroadcastall2 = 0;
auto treduceall2 = 0;
auto totherall2 = 0;

// function to print combinations that contain
// one element from each of the given arrays
void FragCube::solveInquirePointQuery(std::vector<std::vector<int> >& arr, int my_rank, int num_dims, std::string output_folder, int num_procs)
{
    // number of arrays
    int n = arr.size();

    // to keep track of next element in each of
    // the n arrays
    int* indices = new int[n];

    // initialize with first element's index
    for (int i = 0; i < n; i++)
        indices[i] = 0;

    while (1) {

    	std::vector<int> q;

        // print current combination
        for (int i = 0; i < n; i++){
        	q.push_back(arr[i][indices[i]]);
        }
        my_queries.push_back(q);

        // find the rightmost array that has more
        // elements left after the current element
        // in that array
        int next = n - 1;
        while (next >= 0 &&
              (indices[next] + 1 >= arr[next].size())){
            next--;
        }

        // no such array is found so no more
        // combinations left
        if (next < 0){
            return;
        }
        // if found move to next element in that
        // array
        indices[next]++;

        // for all arrays to the right of this
        // array current index again points to
        // first element
        for (int i = next + 1; i < n; i++)
            indices[i] = 0;
    }
}

std::vector<int> intersection2(std::vector<int> const& left_vector, std::vector<int> const& right_vector) {
    auto left = left_vector.begin();
    auto left_end = left_vector.end();
    auto right = right_vector.begin();
    auto right_end = right_vector.end();

    assert(std::is_sorted(left, left_end));
    assert(std::is_sorted(right, right_end));

    std::vector<int> result;

    while (left != left_end && right != right_end) {
        if (*left == *right) {
            result.push_back(*left);
            ++left;
            ++right;
            continue;
        }

        if (*left < *right) {
            ++left;
            continue;
        }

        assert(*left > *right);
        ++right;
    }

    return result;
}

std::vector <int> getIntersectionC2(std::vector < std::vector <int> > &sets){
	std::vector <int> last_intersection_tids = sets[0];
	std::vector<int> curr_intersection_tids;

	for (std::size_t i = 1; i < sets.size(); ++i) {
		curr_intersection_tids = intersection2(last_intersection_tids, sets[i]);
		std::swap(last_intersection_tids, curr_intersection_tids);
		curr_intersection_tids.clear();
	}

	return last_intersection_tids;
}

void FragCube::solvePointQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){
	int result;
	float sum = 0, result_sum;
	std::vector<int> tmp_result;
	std::vector<int> buffer;
	std::vector<int> dimsubset;
	std::vector<int> dimvalues;
	std::vector<
			std::pair<std::vector<int>,
					std::map<std::vector<int>, std::vector<int>>>>::const_iterator it;
	std::map<std::vector<int>, std::vector<int>>::const_iterator it2;
	std::vector<std::vector<int>> sets;

	std::chrono::steady_clock::time_point begin_compute; //Começou a computar consulta

	std::chrono::steady_clock::time_point end_compute; //Terminou de computar consulta

	std::chrono::steady_clock::time_point begin_reduce; //Começou o reduce

	std::chrono::steady_clock::time_point end_reduce; //Terminou o reduce

	begin_compute = std::chrono::steady_clock::now();
	if (std::find(q.begin(), q.end(), -2) == q.end()) { //no inquires
		for (unsigned int i = 0; i < q.size(); i++) {
			if (q[i] != -1) {
				dimsubset.push_back(i);
				dimvalues.push_back(q[i]);
			}
		}

		it =
				std::find_if(cuboids.begin(), cuboids.end(),
						[&dimsubset](const std::pair<std::vector<int>,std::map<std::vector<int>, std::vector<int>>>& element)
						{	return element.first == dimsubset;});

		if (it != cuboids.end()) {
			it2 = it->second.find(dimvalues);
			if (it2 != it->second.end()) {
				sets.push_back(it2->second);
			} else {
				sets.push_back(std::vector<int>());
			}
		} else {
			for (unsigned int i = 0; i < q.size(); i++) {
				if (q[i] != -1) {
					std::vector<int> tmpdim, tmpdim_value;
					tmpdim.push_back(i);
					tmpdim_value.push_back(q[i]);
					it =
							std::find_if(cuboids.begin(), cuboids.end(),
									[&tmpdim](const std::pair<std::vector<int>,std::map<std::vector<int>, std::vector<int>>>& element)
									{	return element.first == tmpdim;});
					it2 = it->second.find(tmpdim_value);
					if (it2 != it->second.end()) {
						sets.push_back(it2->second);
					} else {
						sets.push_back(std::vector<int>());
					}
				}
			}

		}
		if (sets.empty()) {
			for (std::unordered_map<int, std::vector<float>>::iterator it =
					imeas.begin(); it != imeas.end(); ++it) {
				tmp_result.push_back(it->first);
			}
		}

		else if (sets.size() == 1) {
			tmp_result.assign(sets.front().begin(), sets.front().end());
		}

		else {

			std::set_intersection(sets[0].begin(), sets[0].end(),
					sets[1].begin(), sets[1].end(),
					std::back_inserter(tmp_result));

			if (sets.size() > 2) {

				for (size_t i = 2; i < sets.size(); ++i) {
					buffer.clear();
					std::set_intersection(tmp_result.begin(), tmp_result.end(),
							sets[i].begin(), sets[i].end(),
							std::back_inserter(buffer));

					swap(tmp_result, buffer);
				}
			}
		}

		int my_result = tmp_result.size();

		begin_reduce = std::chrono::steady_clock::now();
		MPI_Reduce(&my_result, &result, 1, MPI_FLOAT, MPI_SUM, 0,
				MPI_COMM_WORLD);
		end_reduce = std::chrono::steady_clock::now();

		end_compute = std::chrono::steady_clock::now();

		if (my_rank == 0) {
			auto t1 = std::chrono::duration_cast<std::chrono::microseconds> (end_compute - begin_compute).count();
			auto t2 = std::chrono::duration_cast<std::chrono::microseconds> (end_reduce - begin_reduce).count();
			if (result > 0) {
				for (auto i : q) {
					if (i == -1) {
						std::cout << '*' << ' ';
					} else {
						std::cout << i << ' ';
					}

				}
				//std::cout << ": " << result << " t=" << t1 << "[µs]" <<
				//		" reduce=" << t2 << "[µs]" <<
				//		" other=" << t1 - t2 << "[µs]" << std::endl;
				std::cout << ": " << result << std::endl;
				treduceall2 += t2;
				totherall2 += (t1 - t2);
			}
		}
	} else {
		std::cout << "Query has at least one inquire" << std::endl;
	}
}

void FragCube::solveInquireQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){

	//Listas de atributos em cada dimensão possíveis para esse inquire
	std::vector<std::vector<int>> lists_of_attributes;

	//Listas de TIDs fazem parte da resposta da consulta
	std::vector<std::vector<int>> lists_of_tids;

	//A intereseção final da lista anterior - apenas os BIDs comuns a todas as listas
	std::set<int> last_intersection_tids;

	std::chrono::steady_clock::time_point begin_broadcast; //Começou o broadcast

	std::chrono::steady_clock::time_point end_broadcast; //Terminou o broadcast

	//Percorra a consulta e, quando não encontrar agregação ou inquire, pegue a lista de TIDs da memória
	for(std::size_t i = 0; i < q.size(); i++){
		if (q[i] != -1 && q[i] != -2)
		{
			lists_of_tids.push_back(iindex[i][q[i]]);
		}
	}

	//Guardamos na variável o valor da interseção de todas as listas
    copy(lists_of_tids[0].begin(),lists_of_tids[0].end(),inserter(last_intersection_tids, last_intersection_tids.end()));
	std::set<int> curr_intersection_tids;

	//https://stackoverflow.com/questions/25505868/the-intersection-of-multiple-sorted-arrays
	for (std::size_t i = 1; i < lists_of_tids.size(); ++i) {
		std::set_intersection(last_intersection_tids.begin(), last_intersection_tids.end(),
			lists_of_tids[i].begin(), lists_of_tids[i].end(),
			std::inserter(curr_intersection_tids, curr_intersection_tids.begin()));
		swap(last_intersection_tids, curr_intersection_tids);
		curr_intersection_tids.clear();
	}

	for(std::size_t i = 0; i < q.size(); i++){
		std::vector<int> attribs;
		if (q[i] == -2)
		{
			attribs.push_back(-1);
			for (auto& it: iindex[i]) {
				//std::cout << it.first << " ";
				//attribs.push_back(it.first);

				std::set<int> curr_intersection_tids, curr_tids;
			    copy(it.second.begin(),it.second.end(),inserter(curr_tids, curr_tids.end()));

				std::set_intersection(curr_tids.begin(), curr_tids.end(),
						last_intersection_tids.begin(), last_intersection_tids.end(),
					std::inserter(curr_intersection_tids, curr_intersection_tids.begin()));

				if(curr_intersection_tids.size() > 0){
					attribs.push_back(it.first);
				}

				curr_intersection_tids.clear();
			}
			sort( attribs.begin(), attribs.end() );
			attribs.erase( unique( attribs.begin(), attribs.end() ), attribs.end() );
			lists_of_attributes.insert(lists_of_attributes.begin() + i, attribs);
		}
		else{
			attribs.push_back(q[i]);
			lists_of_attributes.insert(lists_of_attributes.begin() + i, attribs);
		}
	}

	/*for (int i = 0; i < lists_of_attributes.size(); i++)
	{
	    for (int j = 0; j < lists_of_attributes[i].size(); j++)
	    {
	        std::cout << lists_of_attributes[i][j] << " ";
	    }
	    std::cout << "\n";
	}*/


	solveInquirePointQuery(lists_of_attributes, my_rank, num_dims, output_folder, num_procs);

	MPI_Barrier(MPI_COMM_WORLD);

	tbroadcastall2 = 0;

	for(int i = 0; i < num_procs; i++){

		int num_queries = my_queries.size();
		std::vector<int> q;
		q.resize(num_dims);

		begin_broadcast = std::chrono::steady_clock::now();
		MPI_Bcast(&num_queries, 1, MPI_INT, i, MPI_COMM_WORLD);
		end_broadcast = std::chrono::steady_clock::now();

		tbroadcastall2 += std::chrono::duration_cast<std::chrono::microseconds> (end_broadcast - begin_broadcast).count();

		for(int j=0; j<num_queries;j++){
			if(my_rank == i){
				q = my_queries[j];
			}

			begin_broadcast = std::chrono::steady_clock::now();
			MPI_Bcast(&q[0], num_dims, MPI_INT, i, MPI_COMM_WORLD);
			end_broadcast = std::chrono::steady_clock::now();

			tbroadcastall2 += std::chrono::duration_cast<std::chrono::microseconds> (end_broadcast - begin_broadcast).count();

			if(my_rank != i){
				for(int k=0; k<my_queries.size(); k++){
					if(my_queries[k] == q){
						my_queries.erase(my_queries.begin() + k);
					}
				}
			}

			solvePointQuery(q, my_rank, num_dims, output_folder, num_procs);
		}
	}

	if(my_rank == 0){
		std::cout << "Time broadcasting = " << tbroadcastall2 << "[µs]" << std::endl;
		std::cout << "Time reducing = " << treduceall2 << "[µs]" << std::endl;
		std::cout << "Time other = " << totherall2 << "[µs]" << std::endl;
	}

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

void FragCube::ComputeCube(std::string cube_table, int num_dims, int num_meas,
                int partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand)
{

        MPI_File data; //Um handler de arquivos MPI para o dataset completo
        MPI_Offset tid_offset; //Uma forma de saber o endereço de um TID numa tupla do dataset
        MPI_Offset dims_offset; //Uma forma de saber o endereço das dimensões numa tupla do dataset
        MPI_Offset meas_offset; //Uma forma de saber o endereço das medidas numa tupla do dataset

        //Abra o arquivo fornecido pelo usuário e guarde o endereço dele no handler
        MPI_File_open(MPI_COMM_WORLD, cube_table.c_str(), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &data);

        int count = 0;
        int increment = partition_size / reading_rate;
        int fragment_size = 1;
        iindex.resize(num_dims);
        read_buffer.resize(increment);

        //Ajusta o tamanho das tuplas no buffer de leitura, que a princípio não se sabe
        //Assim é possível ler tuplas com a quantidade correta de dimensões e medidas
        for (TupleType &i : read_buffer)
        {
                i.dims.resize(num_dims);
                i.meas.resize(num_meas);
        }

        std::vector<int> dimsubset;
        std::vector<int> idims(num_dims);
        std::vector<std::vector<int>> allDimVecs;
        std::vector<int> DimVec;
        std::vector<int> CuboidCell;
        std::vector<int> tmp_result;
        std::vector<int> buffer;
        std::vector<std::vector<int>> sets;
        std::iota(std::begin(idims), std::end(idims), 0); // Fill with 0, 1, ..., num_dims.
        std::unordered_map<int, std::vector<int>>::const_iterator iter;

        while (count < partition_size)
        {
                if ((count + increment) > partition_size)
                {
                        increment = partition_size - count;
                        read_buffer.resize(increment);
                }

                for (int a = 0; a < increment; ++a)
                {
                        tid_offset =
                                        ((my_rank * partition_size) + (a + count))
                                                        * (sizeof(int)
                                                                        + (num_dims
                                                                                        * sizeof(int))
                                                                        + (num_meas
                                                                                        * sizeof(float)));
                        dims_offset = tid_offset + sizeof(int);
                        meas_offset = dims_offset
                                        + ((num_dims) * sizeof(int));
                        MPI_File_read_at(data, tid_offset, &read_buffer[a], 1,
                        MPI_INT,
                        MPI_STATUS_IGNORE);
                        MPI_File_read_at(data, dims_offset,
                                        &read_buffer[a].dims[0], 1,
                                        MPI_TUPLE_DIMS, MPI_STATUS_IGNORE);
                        MPI_File_read_at(data, meas_offset,
                                        &read_buffer[a].meas[0], 1,
                                        MPI_TUPLE_MEAS,
                                        MPI_STATUS_IGNORE);

                }

                for (TupleType &t : read_buffer)
                {
                        bool use_tuple = true;

                        if (on_demand)
                        {
                                int useful_in_queries = 0;
                                for (std::vector<int> &q : queries)
                                {
                                        if (useful_in_queries == 0)
                                        {
                                                int useful_values = 0;
                                                for (int x = 0; x < num_dims;
                                                                x++)
                                                {
                                                        if ((t.dims[x] == q[x])
                                                                        || q[x]
                                                                                        == -1
                                                                        || q[x]
                                                                                        == -2)
                                                        {
                                                                useful_values++;
                                                        }
                                                }
                                                if (useful_values == num_dims)
                                                {
                                                        useful_in_queries++;
                                                }
                                        }
                                }
                                if (useful_in_queries == 0)
                                {
                                        use_tuple = false;
                                }
                        }

                        if (use_tuple)
                        {
                                for (int j = 0; j < num_dims; ++j)
                                {
                                        iindex[j][t.dims[j]].push_back(t.tid);
                                }

                                for (int k = 0; k < num_meas; ++k)
                                {
                                        imeas[t.tid].push_back(t.meas[k]);
                                }
                        }

                }

                increment = partition_size / reading_rate;
                count += increment;

        }

        int dims_proc = 0;

        while (dims_proc < num_dims)
        {
                for (int s = 1; s < (1 << fragment_size); ++s) // iterate through all non-null sets
                {
                        for (int e = 0; e < fragment_size; ++e) // for each set element
                        {
                                if (s & (1 << e)) // test for membership of set
                                {
                                        dimsubset.push_back(
                                                        idims[e + dims_proc]);
                                }
                        }

                        //Hold the cells of this cuboid
                        std::map<std::vector<int>, std::vector<int>> cuboid;

                        for (auto num2 : dimsubset)
                        {
                                for (auto &kv : iindex[num2])
                                {
                                        //std::cout << kv.first << ' ';
                                        DimVec.push_back(kv.first);
                                }
                                allDimVecs.push_back(DimVec);
                                DimVec.clear();
                        }

                        size_t max = 1;

                        for (auto const &v : allDimVecs)
                                max *= v.size();

                        for (size_t i = 0; i < max; i++)
                        {
                                auto temp = i;
                                for (auto const &vec : allDimVecs)
                                {
                                        auto index = temp % vec.size();
                                        temp /= vec.size();
                                        CuboidCell.push_back(vec[index]);
                                }

                                for (unsigned int z = 0; z < dimsubset.size();
                                                z++)
                                {
                                        if (CuboidCell[z] != -1)
                                        {
                                                iter =
                                                                iindex[dimsubset[z]].find(
                                                                                CuboidCell[z]);
                                                if (iter
                                                                != iindex[dimsubset[z]].end())
                                                {
                                                        sets.push_back(
                                                                        iter->second);
                                                }
                                                else
                                                {
                                                        sets.push_back(
                                                                        std::vector<
                                                                                        int>());
                                                }
                                        }
                                }

                                if (sets.empty())
                                {
                                        for (std::unordered_map<int,
                                                        std::vector<float>>::iterator it =
                                                        imeas.begin();
                                                        it != imeas.end();
                                                        ++it)
                                        {
                                                tmp_result.push_back(
                                                                it->first);
                                        }
                                }

                                else if (sets.size() == 1)
                                {
                                        tmp_result.assign(
                                                        sets.front().begin(),
                                                        sets.front().end());
                                }

                                else
                                {

                                        std::set_intersection(sets[0].begin(),
                                                        sets[0].end(),
                                                        sets[1].begin(),
                                                        sets[1].end(),
                                                        std::back_inserter(
                                                                        tmp_result));

                                        if (sets.size() > 2)
                                        {

                                                for (size_t i = 2;
                                                                i
                                                                                < sets.size();
                                                                ++i)
                                                {
                                                        buffer.clear();
                                                        std::set_intersection(
                                                                        tmp_result.begin(),
                                                                        tmp_result.end(),
                                                                        sets[i].begin(),
                                                                        sets[i].end(),
                                                                        std::back_inserter(
                                                                                        buffer));

                                                        swap(tmp_result,
                                                                        buffer);
                                                }
                                        }
                                }

                                std::vector<int> s(tmp_result.begin(),
                                                tmp_result.end());

                                //insert cuboid into map of cuboids
                                if(!s.empty()){
                                cuboid[CuboidCell] = s;
                                }

                                tmp_result.clear();
                                buffer.clear();
                                sets.clear();
                                CuboidCell.clear();
                        }

                        //cuboid inserted into general list of cubes
                        cuboids.push_back(std::make_pair(dimsubset, cuboid));
                        allDimVecs.clear();
                        dimsubset.clear();
                }
                dims_proc += fragment_size;
                if (num_dims - dims_proc < fragment_size)
                {
                        fragment_size = num_dims - dims_proc;
                }
        }

        read_buffer.clear();

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_close(&data);

}

void FragCube::QueryCube(std::vector<std::vector<int>> queries, int my_rank, int num_dims, std::string output_folder,  int num_procs){

	//Para cada consulta faça
    for (std::vector<int> &q : queries)
    {
        if (std::find(q.begin(), q.end(), -2) == q.end())
        { //no inquires
        	solvePointQuery(q, my_rank, num_dims, output_folder, num_procs);
        } else {//pelo menos um inquire
        	solveInquireQuery(q, my_rank, num_dims, output_folder, num_procs);
        }
    }
}
