#include "BlockCube.h"
#include <algorithm>

#include <chrono>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

auto tbroadcastall = 0;
auto treduceall = 0;
auto totherall = 0;

// function to print combinations that contain
// one element from each of the given arrays
void BlockCube::solveInquirePointQuery(std::vector<std::vector<int> >& arr, int my_rank, int num_dims, std::string output_folder, int num_procs)
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
              (static_cast<std::vector<int>::size_type>(indices[next] + 1) >= arr[next].size())){
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

std::vector<int> intersection(std::vector<int> const& left_vector, std::vector<int> const& right_vector) {
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

std::vector <int> getIntersectionC(std::vector < std::vector <int> > &sets){
	std::vector <int> last_intersection_tids = sets[0];
	std::vector<int> curr_intersection_tids;

	for (std::size_t i = 1; i < sets.size(); ++i) {
		curr_intersection_tids = intersection(last_intersection_tids, sets[i]);
		std::swap(last_intersection_tids, curr_intersection_tids);
		curr_intersection_tids.clear();
	}

	return last_intersection_tids;
}

std::vector <int> getIntersectionA(std::vector < std::vector <int> > &sets){
	std::vector <int> last_intersection_tids = sets[0];
	std::vector<int> curr_intersection_tids;

	for (std::size_t i = 1; i < sets.size(); ++i) {
		std::set_intersection(last_intersection_tids.begin(), last_intersection_tids.end(),
				sets[i].begin(), sets[i].end(),
			std::back_inserter(curr_intersection_tids));
		std::swap(last_intersection_tids, curr_intersection_tids);
		curr_intersection_tids.clear();
	}

	return last_intersection_tids;
}

std::vector <int> getIntersectionB(std::vector < std::vector <int> > &sets)
{
	std::vector <int> result;  // To store the reaultant set
int smallSetInd = 0;  // Initialize index of smallest set
int minSize = sets[0].size(); // Initialize size of smallest set

// sort all the sets, and also find the smallest set
for (int i = 1 ; static_cast<std::vector<int>::size_type>(i) < sets.size() ; i++)
{
    // sort this set
    //sort(sets[i].begin(), sets[i].end());

    // update minSize, if needed
    if (static_cast<std::vector<int>::size_type>(minSize) > sets[i].size())
    {
        minSize = sets[i].size();
        smallSetInd = i;
    }
}

std::map<int,int> elementsMap;

// Add all the elements of smallest set to a map, if already present,
// update the frequency
for (int i = 0; static_cast<std::vector<int>::size_type>(i) < sets[smallSetInd].size(); i++)
{
    if (elementsMap.find( sets[smallSetInd][i] ) == elementsMap.end())
        elementsMap[ sets[smallSetInd][i] ] = 1;
    else
        elementsMap[ sets[smallSetInd][i] ]++;
}

// iterate through the map elements to see if they are present in
// remaining sets
std::map<int,int>::iterator it;
for (it = elementsMap.begin(); it != elementsMap.end(); ++it)
{
    int elem = it->first;
    int freq = it->second;

    bool bFound = true;

    // Iterate through all sets
    for (int j = 0 ; static_cast<std::vector<int>::size_type>(j) < sets.size() ; j++)
    {
        // If this set is not the smallest set, then do binary search in it
        if (j != smallSetInd)
        {
            // If the element is found in this set, then find its frequency
            if (binary_search( sets[j].begin(), sets[j].end(), elem ))
            {
               int lInd = lower_bound(sets[j].begin(), sets[j].end(), elem)
                                                        - sets[j].begin();
               int rInd = upper_bound(sets[j].begin(), sets[j].end(), elem)
                                                        - sets[j].begin();

               // Update the minimum frequency, if needed
               if ((rInd - lInd) < freq)
                   freq = rInd - lInd;
            }
            // If the element is not present in any set, then no need
            // to proceed for this element.
            else
            {
                bFound = false;
                break;
            }
        }
    }

    // If element was found in all sets, then add it to result 'freq' times
    if (bFound)
    {
        for (int k = 0; k < freq; k++)
            result.push_back(elem);
    }
}
return result;
}

void BlockCube::solvePointQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){
	//Listas de BIDs que são úteis para responder a consulta
	std::vector<std::set<int>> lists_of_bids;

	//A intereseção final da lista anterior - apenas os BIDs comuns a todas as listas
	std::set<int> last_intersection_bids;

	//Listas de TIDs fazem parte da resposta da consulta
	std::vector<std::vector<int>> lists_of_tids;

	//A intereseção final da lista anterior - apenas os TIDs comuns a todas as listas
	std::vector<int> last_intersection_tids;

	std::vector<std::vector<int>> lists_of_attributes;

	int my_result;

	int result;

	std::chrono::steady_clock::time_point begin_compute; //Começou a computar consulta

	std::chrono::steady_clock::time_point end_compute; //Terminou de computar consulta

	std::chrono::steady_clock::time_point begin_reduce; //Começou o reduce

	std::chrono::steady_clock::time_point end_reduce; //Terminou o reduce

	begin_compute = std::chrono::steady_clock::now();
	//Percorra a consulta e, quando não encontrar agregação, pegue a lista de BIDs da memória
	for(std::size_t i = 0; i < q.size(); i++){
		if (q[i] != -1)
		{
			lists_of_bids.push_back(bblocRAM[i][q[i]]);
		}
	}

	//Guardamos na variável o valor da interseção de todas as listas
	last_intersection_bids = lists_of_bids[0];
	std::set<int> curr_intersection_bids;

	//https://stackoverflow.com/questions/25505868/the-intersection-of-multiple-sorted-arrays
	for (std::size_t i = 1; i < lists_of_bids.size(); ++i) {
		std::set_intersection(last_intersection_bids.begin(), last_intersection_bids.end(),
			lists_of_bids[i].begin(), lists_of_bids[i].end(),
			std::inserter(curr_intersection_bids, curr_intersection_bids.begin()));
		swap(last_intersection_bids, curr_intersection_bids);
		curr_intersection_bids.clear();
	}

	if(last_intersection_bids.size() > 0){
		for(int i = 0 ; i < num_dims ; i++) bbloc[i].resize(*last_intersection_bids.rbegin()+1);

		for(auto bid : last_intersection_bids){
			for(int i = 0; i < num_dims; i++){
				if (q[i] != -1 && bbloc[i][bid].empty()){
					std::string filename = output_folder + "/" + std::to_string(my_rank) + "/" + std::to_string(i) + "/" + std::to_string(bid);
					std::ifstream ifs(filename.c_str(), std::ifstream::binary);
					boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);
					ia & bbloc[i][bid];
				}
			}
		}

		for(std::size_t i = 0; i < q.size(); i++){
			std::vector<int> tids;
			if (q[i] != -1)
			{
				for(auto bid : last_intersection_bids){
					tids.insert(std::end(tids), std::begin(bbloc[i][bid][q[i]]), std::end(bbloc[i][bid][q[i]]));
				}
				lists_of_tids.push_back(tids);
			}
		}

		//Ordena as listas de tid por tamanho
		std::sort(lists_of_tids.begin(), lists_of_tids.end(), [](const std::vector<int> & a, const std::vector<int> & b){ return a.size() < b.size(); });

		last_intersection_tids = getIntersectionC(lists_of_tids);

	}

	my_result = last_intersection_tids.size();

	begin_reduce = std::chrono::steady_clock::now();
	MPI_Reduce(&my_result, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	end_reduce = std::chrono::steady_clock::now();

	end_compute = std::chrono::steady_clock::now();

	if (my_rank == 0)
	{
			auto t1 = std::chrono::duration_cast<std::chrono::microseconds> (end_compute - begin_compute).count();
			auto t2 = std::chrono::duration_cast<std::chrono::microseconds> (end_reduce - begin_reduce).count();
			if (result > 0)
			{
					for (auto i : q)
					{
							if (i == -1)
							{
									std::cout << '*'
													<< ' ';
							}
							else
							{
									std::cout << i << ' ';
							}

					}
					//std::cout << ": " << result << " t=" << t1 << "[µs]" <<
					//		" reduce=" << t2 << "[µs]" <<
					//		" other=" << t1 - t2 << "[µs]" << std::endl;
					std::cout << ": " << result << std::endl;
					treduceall += t2;
					totherall += (t1 - t2);
			}
	}
}

void BlockCube::solveInquireQuery(std::vector<int> q, int my_rank, int num_dims, std::string output_folder, int num_procs){
	//Listas de BIDs que são úteis para responder a consulta
	std::vector<std::set<int>> lists_of_bids;

	//A intereseção final da lista anterior - apenas os BIDs comuns a todas as listas
	std::set<int> last_intersection_bids;

	//Listas de TIDs fazem parte da resposta da consulta
	std::vector<std::vector<int>> lists_of_tids;

	//A intereseção final da lista anterior - apenas os TIDs comuns a todas as listas
	std::vector<int> last_intersection_tids;

	//Listas de atributos em cada dimensão possíveis para esse inquire
	std::vector<std::vector<int>> lists_of_attributes;

	std::chrono::steady_clock::time_point begin_broadcast; //Começou o broadcast

	std::chrono::steady_clock::time_point end_broadcast; //Terminou o broadcast

	//Percorra a consulta e, quando não encontrar agregação ou inquire, pegue a lista de BIDs da memória
	for(std::size_t i = 0; i < q.size(); i++){
		if (q[i] != -1 && q[i] != -2)
		{
			lists_of_bids.push_back(bblocRAM[i][q[i]]);
		}
	}

	//Guardamos na variável o valor da interseção de todas as listas
	last_intersection_bids = lists_of_bids[0];
	std::set<int> curr_intersection_bids;

	//https://stackoverflow.com/questions/25505868/the-intersection-of-multiple-sorted-arrays
	for (std::size_t i = 1; i < lists_of_bids.size(); ++i) {
		std::set_intersection(last_intersection_bids.begin(), last_intersection_bids.end(),
			lists_of_bids[i].begin(), lists_of_bids[i].end(),
			std::inserter(curr_intersection_bids, curr_intersection_bids.begin()));
		swap(last_intersection_bids, curr_intersection_bids);
		curr_intersection_bids.clear();
	}


	for(int i = 0 ; i < num_dims ; i++) bbloc[i].resize(*last_intersection_bids.rbegin()+1);

	for(auto bid : last_intersection_bids){
		for(int i = 0; i < num_dims; i++){
			if(q[i] == -2){
				std::string filename = output_folder + "/" + std::to_string(my_rank) + "/" + std::to_string(i) + "/" + std::to_string(bid);
				std::ifstream ifs(filename.c_str(), std::ifstream::binary);
				boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);
				ia & bbloc[i][bid];
			}
		}
	}

	for(std::size_t i = 0; i < q.size(); i++){
		std::vector<int> attribs;
		if (q[i] == -2)
		{
			attribs.push_back(-1);
			for(auto bid : last_intersection_bids){
				//std::cout << "ATRIBUTOS DA DIMENSÃO " << i << " NO BID " << bid << std::endl;
				for (auto& it: bbloc[i][bid]) {
					//std::cout << it.first << " ";
					attribs.push_back(it.first);
				}
				//std::cout << std::endl;
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

	tbroadcastall = 0;

	for(int i = 0; i < num_procs; i++){

		int num_queries = my_queries.size();
		std::vector<int> q;
		q.resize(num_dims);

		begin_broadcast = std::chrono::steady_clock::now();
		MPI_Bcast(&num_queries, 1, MPI_INT, i, MPI_COMM_WORLD);
		end_broadcast = std::chrono::steady_clock::now();

		tbroadcastall += std::chrono::duration_cast<std::chrono::microseconds> (end_broadcast - begin_broadcast).count();

		for(int j=0; j<num_queries;j++){
			if(my_rank == i){
				q = my_queries[j];
			}

			begin_broadcast = std::chrono::steady_clock::now();
			MPI_Bcast(&q[0], num_dims, MPI_INT, i, MPI_COMM_WORLD);
			end_broadcast = std::chrono::steady_clock::now();

			tbroadcastall += std::chrono::duration_cast<std::chrono::microseconds> (end_broadcast - begin_broadcast).count();

			if(my_rank != i){
				for(int k=0; static_cast<std::vector<int>::size_type>(k)<my_queries.size(); k++){
					if(my_queries[k] == q){
						my_queries.erase(my_queries.begin() + k);
					}
				}
			}

			solvePointQuery(q, my_rank, num_dims, output_folder, num_procs);
		}
	}


	if(my_rank == 0){
		std::cout << "Time broadcasting = " << tbroadcastall << "[µs]" << std::endl;
		std::cout << "Time reducing = " << treduceall << "[µs]" << std::endl;
		std::cout << "Time other = " << totherall << "[µs]" << std::endl;
	}

}

BlockCube::BlockCube()
{
}

BlockCube::BlockCube(const BlockCube&)
{
}

BlockCube::~BlockCube()
{
}

BlockCube&
BlockCube::operator=(const BlockCube&)
{
        return *this;
}

void BlockCube::ComputeCube(std::string cube_table, int num_dims,
                int num_meas, int partition_size, int reading_rate, int tbloc, std::string output_folder, int my_rank,
                std::vector<std::vector<int>> queries, bool on_demand)
{

        MPI_File data; //Um handler de arquivos MPI para o dataset completo
        MPI_Offset tid_offset; //Uma forma de saber o endereço de um TID numa tupla do dataset
        MPI_Offset dims_offset; //Uma forma de saber o endereço das dimensões numa tupla do dataset
        MPI_Offset meas_offset; //Uma forma de saber o endereço das medidas numa tupla do dataset

        //Abra o arquivo fornecido pelo usuário e guarde o endereço dele no handler
        MPI_File_open(MPI_COMM_WORLD, cube_table.c_str(), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &data);

        int tuples_read = 0; //Conta quantas tuplas foram lidas
        int buffer_read = 0; //Conta quantas vezes o buffer de leitura foi preenchido
        int read_increment = partition_size / reading_rate; //Lê a partição com base na taxa de ingestão
        int bid = 0; //BID inicial
        int bid_tuples = 0; //Quantidade de tuplas no bloco identificado por BID

        //Na definição da classe essas estruturas não tem tamanho definido ainda
        //Aqui podemos redimensionar elas de acordo com o dataset
        bbloc.resize(num_dims);
        bblocRAM.resize(num_dims);
        read_buffer.resize(read_increment);
        for(int i = 0 ; i < num_dims ; i++) bbloc[i].resize(1); //Aloca espaço para o BID inicial em todas as dimensões

        //Ajusta o tamanho das tuplas no buffer de leitura, que a princípio não se sabe
        //Assim é possível ler tuplas com a quantidade correta de dimensões e medidas
        for (TupleType &tuple : read_buffer)
        {
        	tuple.dims.resize(num_dims);
        	tuple.meas.resize(num_meas);
        }

        //Creates a directory in the output folder for the current process (bcubing data will be stored inside)
        boost::filesystem::create_directories(output_folder + "/" + std::to_string(my_rank));

        //Creates directories for each dimension of the dataset on the process folder, initially empty (will hold block data)
        for (int i = 0; i < num_dims; i++)
        {
                boost::filesystem::create_directories(output_folder + "/" + std::to_string(my_rank) + "/" + std::to_string(i));
        }

        //Enquanto não ler todas as tuplas dessa partição
        while (tuples_read < partition_size)
        {
        		//Verifica se na próxima leitura que vai fazer não irá passar do tamanho máximo da partição
        		//Se necessário, reduz o incremento e o buffer de leitura
                if ((tuples_read + read_increment) > partition_size)
                {
                        read_increment = partition_size - tuples_read;
                        read_buffer.resize(read_increment);
                }

                //Lê as tuplas do incremento atual e salva no buffer de leitura
                for (int next_tuple = 0; next_tuple < read_increment; ++next_tuple)
                {
                		//A posição do próximo TID no arquivo binário
                        tid_offset =
                                        ((my_rank * partition_size) + (next_tuple + tuples_read))
                                                        * (sizeof(int)
                                                                        + (num_dims
                                                                                        * sizeof(int))
                                                                        + (num_meas
                                                                                        * sizeof(float)));
                        //A posição das próximas dimensões no arquivo binário
                        dims_offset = tid_offset + sizeof(int);

                        //A posição das próximas medidas no arquivo binário
                        meas_offset = dims_offset
                                        + ((num_dims) * sizeof(int));

                        //Leia um inteiro e guarde no buffer
                        MPI_File_read_at(data, tid_offset, &read_buffer[next_tuple], 1,
                        MPI_INT,
                        MPI_STATUS_IGNORE);

                        //Leia todas as dimensões e guarde no buffer
                        MPI_File_read_at(data, dims_offset,
                                        &read_buffer[next_tuple].dims[0], 1,
                                        MPI_TUPLE_DIMS, MPI_STATUS_IGNORE);

                        //Leia todas as medidas e guarde no buffer
                        MPI_File_read_at(data, meas_offset,
                                        &read_buffer[next_tuple].meas[0], 1,
                                        MPI_TUPLE_MEAS,
                                        MPI_STATUS_IGNORE);

                }

                //Nesse ponto o buffer de leitura foi preenchido

                //Para cada uma das tuplas do buffer
                for (TupleType &t : read_buffer)
                {
                		//A princípio considere que a tupla é "útil"
                        bool use_tuple = true;

                        //Se cubo sob demanda, verifique se a tupla é realmente "útil"
                        if (on_demand)
                        {
                        		//Quantidade de consultas que envolvem a tupla
                                int useful_in_queries = 0;

                                //Para cada consulta a ser executada
                                for (std::vector<int> &q : queries)
                                {
                                        if (useful_in_queries == 0)
                                        {
                                        		//Valores úteis na tupla para a consulta
                                                int useful_values = 0;

                                                //Veja se cada valor da tupla bate com a consulta
                                                //Também valida se consulta tem agregação ou inquire
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

                                                //Se todos os valores atenderem a consulta, então a tupla é "útil"
                                                if (useful_values == num_dims)
                                                {
                                                        useful_in_queries++;
                                                }
                                        }
                                }

                                //Se não for útil para a consulta, marque para não usar
                                if(useful_in_queries == 0){
                                        use_tuple = false;
                                }
                        }

                        //Tupla será inserida no cubo de dados
                        if (use_tuple)
                        {
                        	//Conta uma nova tupla no bloco atual
                        	bid_tuples++;

                        	//Se essa nova tupla significar que irá exceder tamanho do bloco
                        	//Crie um novo BID e coloque a tupla nele
                        	if(bid_tuples > tbloc){
                                for (int i = 0; i < num_dims; i++)
                                {
                                        std::string filename = output_folder + "/" + std::to_string(my_rank) + "/" + std::to_string(i) + "/" + std::to_string(bid);
                                        std::ofstream ofs(filename.c_str(), std::ofstream::binary);
                                        boost::archive::binary_oarchive oa(ofs, boost::archive::no_header);
                                        oa & bbloc[i][0];

                                }
                                bbloc.clear();
                                bbloc.resize(num_dims);
                                for(int i = 0 ; i < num_dims ; i++) bbloc[i].resize(1); //Aloca espaço para o BID inicial em todas as dimensões

                        		bid++;
                        		bid_tuples = 1;//Novo bloco terá uma tupla
                        	}

							for (int j = 0; j < num_dims; ++j)
							{
									bbloc[j][0][t.dims[j]].push_back(t.tid);
									bblocRAM[j][t.dims[j]].insert(bid);
							}

							for (int k = 0; k < num_meas; ++k)
							{
									bmeas[t.tid].push_back(t.meas[k]);
							}

                        }

                }


                for (int i = 0; i < num_dims; i++)
                {
                        std::string filename = output_folder + "/" + std::to_string(my_rank) + "/" + std::to_string(i) + "/" + std::to_string(bid);
                        std::ofstream ofs(filename.c_str(), std::ofstream::binary);
                        boost::archive::binary_oarchive oa(ofs, boost::archive::no_header);
                        oa & bbloc[i][0];

                }

                read_increment = partition_size / reading_rate;
                tuples_read += read_increment;
                buffer_read++;
        }

        bbloc.clear();
        bbloc.resize(num_dims);
        read_buffer.clear();

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_close(&data);
}


void BlockCube::QueryCube(std::vector<std::vector<int>> queries, int my_rank, int num_dims, std::string output_folder,  int num_procs){

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
