#include "BlockCube.h"
#include <algorithm>
#include <numeric>
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

/**
 * Computa um cubo de dados - essa é uma função fixa do framework
 *
 * @param cube_table o nome de arquivo para o dataset
 * @param num_dims número de dimensões do dataset
 * @param num_meas número de medidas do dataset
 * @param tuple_partition_size o tamanho da partição (por tuplas)
 * @param dim_partition_size o tamanho da partição (por dimensões)
 * @param reading_rate a taxa de leitura do disco usada pelos processos
 * @param tbloc variável específica do algoritmo bCubing
 * @param output_folder arquivo de saída, se o cubo for escrito em disco
 * @param my_rank o rank MPI do processo atual
 * @param queries as consultas que serão executadas
 * @param on_demand booleano que determina se a computação será feita "on-demand"eiros.
 *
 */
void BlockCube::ComputeCube(std::string cube_table, int num_dims,
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

        ///////////////////// ESPECÍFICO DO ALGORITMO BCUBING /////////////////////

        //ID do primeiro bloco, blocos subsequentes terão IDs incrementais
        int bid = 0;

        //Quantidade de tuplas no bloco identificado por BID, deve ser igual a tbloc
        int bid_tuples_count = 0;

        //Sabendo o tamanho da partição e o tamanho de bloco conseguiremos saber quantos blocos serão criados
        int num_bids = (tuple_partition_size % tbloc == 0) ? tuple_partition_size / tbloc : tuple_partition_size / tbloc + 1;

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBloc
        bbloc.resize(num_dims);

        //Agora que sabemos a quantidade de dimensões, redimensiona a estrutura bCubingBlocRAM
        bblocRAM.resize(num_dims);

        //Agora que sabemos a quantidade de blocos, redimensiona a estrutura bCubingMeasures
        bmeas.resize(num_bids);

        //Calculamos um bloco por vez apenas, em todas as dimensões, assim utiliza menos memória
        //Depois de processar o bloco será salvo em disco imediatamente
        for (auto & bids_per_dimension : bbloc) {
        	bids_per_dimension.resize(1);
        }

        //Se o diretório de saída estiver vazio, significa que foi recém criado (cubo ainda não computado)
        //Nesse caso procede à computação do cubo normalmente
        if(boost::filesystem::is_empty(output_folder)){

        	//Cada processo MPI tem um diretóri próprio para armazenar dados do cubo
        	std::string process_directory = output_folder + "/" + std::to_string(my_rank);

            //Cria o diretório associado ao processo
            boost::filesystem::create_directory(process_directory);

            //No diretório do processo cria diretórios (inicialmente vazios) para cada uma das dimensões
            //Os blocos usados pelo algoritmo bcubing serão salvos nesses diretórios
            for (int dim_number = 0; dim_number < num_dims; dim_number++)
            {
                    boost::filesystem::create_directory(process_directory + "/" + std::to_string(dim_number));
            }

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

                            //Tupla será inserida no cubo de dados
                            if (use_tuple)
                            {
                            	//Conta uma nova tupla no bloco atual
                            	bid_tuples_count++;

                            	//Se essa nova tupla significar que irá exceder tamanho do bloco
                            	//Crie um novo BID e coloque a tupla nele
                            	if(bid_tuples_count > tbloc){

                            		//Primeiro salva os dados do BID já processado, que estava em memória, para disco
                                    for (int dim_number = 0; dim_number < num_dims; dim_number++)
                                    {
                                    		//SALVA O BLOCO EM DISCO

                                    		//Nome do arquivo onde o BID será salvo - diretório do processo, num diretório específico da dimensão
                                            std::string bid_filename = process_directory + "/" + std::to_string(dim_number) + "/" + std::to_string(bid);

                                            //Indica que o arquivo de saída será um stream de dados binários
                                            std::ofstream ofb(bid_filename.c_str(), std::ofstream::binary);

                                            //Serialização é feita pela biblioteca Boost
                                            boost::archive::binary_oarchive ob(ofb, boost::archive::no_header);

                                            //Salva os dados do BID recém criado em disco
                                            ob & bbloc[dim_number][0];

                                            //SALVA AS MEDIDAS DO BLOCO EM DISCO

                                            //Nome do arquivo onde as medidas do BID serão salvas - diretório do processo, num diretório específico da dimensão
                                            std::string bid_meas_filename = process_directory + "/" + std::to_string(dim_number) + "/" + std::to_string(bid) + ".m";

                                            //Indica que o arquivo de saída será um stream de dados binários
                                            std::ofstream ofm(bid_meas_filename.c_str(), std::ofstream::binary);

                                            //Serialização é feita pela biblioteca Boost
                                            boost::archive::binary_oarchive om(ofm, boost::archive::no_header);

                                            //Salva os dados de medidas do BID recém criado em disco
                                            om & bmeas[0];

                                    }

                                    //Limpa a memória da variável bbloc
                                    bbloc.clear();

                                    //Redimensiona novamente com base na quantidade de dimensões
                                    bbloc.resize(num_dims);

                                    //Limpa a memória da variável bmeas
                                    bmeas.clear();

                                    //Redimensiona novamente com base na quantidade de blocos
                                    bmeas.resize(num_bids);

                                    //Calculamos um bloco por vez apenas, em todas as dimensões, assim utiliza menos memória
                                    //Depois de processar o bloco será salvo em disco imediatamente
                                    for (auto & bids_per_dimension : bbloc) {
                                    	bids_per_dimension.resize(1);
                                    }

                                    //Identificador do novo bloco
                            		bid++;

                            		//Novo bloco já começa com uma tupla, que não foi possível inserir no bloco anterior
                            		bid_tuples_count = 1;
                            	}


                            	//Atualiza as entradas das variáveis bbloc e bblocRAM
                            	//Inclui nessas estruturas os dados da tupla
    							for (int dim_number = 0; dim_number < num_dims; ++dim_number)
    							{
    									//Em cada dimensão, no único bloco existente na posição zero (processamos um bloco por vez),
    									//insere ou atualiza uma entrada com o valor de atributo da dimensão e um novo TID associado
    									bbloc[dim_number][0][tuple.dims[dim_number]].push_back(tuple.tid);

    									//Em cada dimensão insere ou atualiza uma entrada com o valor de atributo da dimensão e um novo BID associado
    									bblocRAM[dim_number][tuple.dims[dim_number]].insert(bid);
    							}

    							//Atualiza as entradas da variável bCubingMeasures
    							//Inclui nessa estrutura os dados da tupla
    							for (int meas_number = 0; meas_number < num_meas; ++meas_number)
    							{
    									//Insere ou atualiza uma entrada com o valor do TID e as medidas associadas
    									//Preenche um bloco por vez na memória e depois salva em disco, por isso aqui sempre usa a primeira posição
    									bmeas[0][tuple.tid].push_back(tuple.meas[meas_number]);
    							}

                            }

                    }

                    //A última escrita de dados da memória pro disco
                    //Ao longo do laço anterior deve ter escrito várias vezes
                    //Aqui é o *último bloco*
                    for (int dim_number = 0; dim_number < num_dims; dim_number++)
                    {
                		//SALVA O BLOCO EM DISCO

                		//Nome do arquivo onde o BID será salvo - diretório do processo, num diretório específico da dimensão
                        std::string bid_filename = process_directory + "/" + std::to_string(dim_number) + "/" + std::to_string(bid);

                        //Indica que o arquivo de saída será um stream de dados binários
                        std::ofstream ofb(bid_filename.c_str(), std::ofstream::binary);

                        //Serialização é feita pela biblioteca Boost
                        boost::archive::binary_oarchive ob(ofb, boost::archive::no_header);

                        //Salva os dados do BID recém criado em disco
                        ob & bbloc[dim_number][0];

                        //SALVA AS MEDIDAS DO BLOCO EM DISCO

                        //Nome do arquivo onde as medidas do BID serão salvas - diretório do processo, num diretório específico da dimensão
                        std::string bid_meas_filename = process_directory + "/" + std::to_string(dim_number) + "/" + std::to_string(bid) + ".m";

                        //Indica que o arquivo de saída será um stream de dados binários
                        std::ofstream ofm(bid_meas_filename.c_str(), std::ofstream::binary);

                        //Serialização é feita pela biblioteca Boost
                        boost::archive::binary_oarchive om(ofm, boost::archive::no_header);

                        //Salva os dados de medidas do BID recém criado em disco
                        om & bmeas[0];

                    }

                    //Pode ser que o incremento de leitura tenha sido alterado acima, então volta para o valor padrão
                    read_increment = tuple_partition_size / reading_rate;

                    //Adiciona a quantidade de tuplas recém lidas ao total de tuplas lidas
                    tuples_read += read_increment;
            }
        }

        //Limpa o espaço em memória da bCubingBloc
        bbloc.clear();

        //Redimensiona pela quantidade de dimensões pois será usado nas consultas
        bbloc.resize(num_dims);

        //Limpa a memória da variável bmeas
        bmeas.clear();

        //Redimensiona novamente com base na quantidade de blocos
        bmeas.resize(num_bids);

        //Limpa o espaço em memória para o buffer de leitura
        read_buffer.clear();

        //Barreira de sincronismo - só continua quando todos os processos chegarem até aqui
        MPI_Barrier(MPI_COMM_WORLD);

        //Todos os processos fecham a conexão com o arquivo de entrada
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
