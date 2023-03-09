#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include "Handler.h"

namespace po = boost::program_options;


/**
 * Converte uma consulta em formato de string para inteiros.
 *
 * @param q_str representa uma consulta ao cubo, composta de caracteres.
 * @return q_int a consulta recebida, convertida em uma sequência de inteiros.
 *
 */
std::vector<int> QueryStoi(std::string q_str) {

	//A consulta, traduzida para inteiros
	std::vector<int> q_int;

	//Armazena um dos argumentos da consulta
	std::string arg;

	//A consulta, em formato de string
	std::stringstream q(q_str);

	//Passa por cada argumento da consulta, separado por espaços
	while (q >> arg) {
		//O símbolo de agregação '*' é traduzido para -1
		if (arg == "*") {
			q_int.push_back(-1);
		//O símbolo de inquire '?' é traduzido para -2
		} else if (arg == "?") {
			q_int.push_back(-2);
		//Converte o argumento para inteiro, quando possível
		} else {
			q_int.push_back(stoi(arg));
		}
	}

	//Retorna a consulta em forma de números inteiros
	return q_int;
}


/**
 * Converte um conjunto de consultas em formato string para inteiros.
 *
 * @param q_vec_str um conjunto de strings com todas as consultas recebidas como parâmetro de entrada.
 * @param q_vec_int um conjunto (espaço em memória) para armazenar as consultas convertidas para inteiros.
 * @param n_dims o número de dimensões das consultas, devem ter o mesmo número de dimensões
 * @return boolean true ou false, se a execução foi bem sucedida ou mal sucedida.
 *
 */
bool QueryConverter(const std::vector<std::string> &q_vec_str,
                std::vector<std::vector<int>> &q_vec_int, int num_dims)
{
		//Passa por cada um dos elementos do conjunto de consultas em formato string
        for (const auto& q_str : q_vec_str)
        {
        		std::vector<int> q_int;

        		//Tenta converter a consulta para inteiros
        		//Só falha se a consulta tiver operadores inválidos ou fora de alcance
                try
                {
                	q_int = QueryStoi(q_str);
                	q_vec_int.push_back(q_int);
                }
                catch (std::invalid_argument& e)
                {
                        std::cout << "ERRO: a consulta " << q_str
                                        << " tem operadores inválidos"
                                        << std::endl;
                        return false;
                }
                catch (std::out_of_range& e)
                {
                        std::cout << "ERRO: a consulta " << q_str
                                        << " tem operadores fora do alcance"
                                        << std::endl;
                        return false;
                }

                //Nesse ponto a consulta foi convertida de string para inteiro
                //No entanto ainda é necessário testar se a quantidade de operadores é diferente
                //da quantidade de dimensões
                if (q_int.size() != static_cast<std::vector<int>::size_type>(num_dims))
                {
                        std::cout << "ERRO: a consulta " << q_str
                                        << " tem uma quantidade de operadores diferente do número de dimensões"
                                        << std::endl;
                        return false;
                }

        }
        return true;
}


/**
 * Calcula o tamanho da partição de um processo.
 *
 * @param num_tuples o número de tuplas no dataset.
 * @param num_procs o número de processos em paralelo sendo executados.
 * @param my_rank o rank do processo atual
 * @param partition_size o tamanho da partição (é calculado dentro dessa função)
 * @param reading_rate a taxa de ingestão de dados do disco, o número de leituras
 * @return boolean true ou false, se a execução foi bem sucedida ou mal sucedida.
 *
 */
bool GetPartitionSize(int num_tuples, int num_dims, int num_procs, int my_rank,
                int &tuple_partition_size, int &dim_partition_size, int reading_rate, std::string partition)
{
		//Usuário escolheu particionar por tuplas
		if(partition == "tuples"){

			//Se houver mais processos em paralelo do que tuplas, é impossível fazer o particionamento
	        if (num_tuples < num_procs)
	        {
	                std::cout
	                                << "ERRO: número de tarefas em paralelo é maior que o número de tuplas!"
	                                << std::endl;
	                return false;
	        }

	        //Se a divisão das tuplas entre os processos é exata, então o tamanho da partição também é exato e igual para todos
	        if (num_tuples % num_procs == 0)
	        {
	        	tuple_partition_size = num_tuples / num_procs;
	        }
	        else{
	            //A divisão não é exata.
	            //Os primeiros processos tem uma tupla a mais
	            if (my_rank < (num_tuples % num_procs))
	            {
	            	tuple_partition_size = (num_tuples / num_procs) + 1;
	            }
	            else
	            {
	            	tuple_partition_size = num_tuples / num_procs;
	            }
	        }

	        //Não é possível ler 1/N de uma quantidade de tuplas menor que N
	        if (tuple_partition_size < reading_rate)
	        {
	                std::cout << "ERRO: a tarefa " << my_rank << " não consegue ler 1/"
	                                << reading_rate << " das " << tuple_partition_size
	                                << " tuplas que ela recebeu!" << std::endl;
	                return false;
	        }

	        //Se escolheu partição por tuplas, isso será ignorado posteriormente
	        dim_partition_size = -1;

	        return true;
		}
		//Usuário escolheu particionar por dimensões
		else if(partition == "dimensions"){

			//Se houver mais processos em paralelo do que tuplas, é impossível fazer o particionamento
	        if(num_dims < num_procs)
	        {
	                std::cout
	                                << "ERRO: número de tarefas em paralelo é maior que o número de dimensões!"
	                                << std::endl;
	                return false;
	        }

	        //Se a divisão das dimensões entre os processos é exata, então o tamanho da partição também é exato e igual para todos
	        if (num_dims % num_procs == 0){
	        	dim_partition_size = num_dims / num_procs;
	        }
	        else{
	            //A divisão não é exata.
	            //Os primeiros processos tem uma dimensão a mais
	            if (my_rank < (num_dims % num_procs))
	            {
	            	dim_partition_size = (num_dims / num_procs) + 1;
	            }
	            else
	            {
	            	dim_partition_size = num_dims / num_procs;
	            }
	        }

	        //Se escolheu partição por dimensões, isso será ignorado posteriormente
	        tuple_partition_size = -1;

	        return true;
		}
		else {
            std::cout
                            << "ERRO: tipo de particionamento não suportado (conhecidos: tuples, dimensions)!"
                            << std::endl;
            return false;
		}
}

/** Lê os dados de entrada e realiza operações adicionaos sobre eles se necessário
*
* @param my_rank o rank MPI do processo atual
* @param num_procs número total de processos sendo executados
* @param num_dims número de dimensões do dataset
* @param num_meas número de medidas do dataset
* @param num_tuples número de tuplas do dataset
* @param tuple_partition_size o tamanho da partição (por tuplas)
* @param dim_partition_size o tamanho da partição (por dimensões)
* @param reading_rate a taxa de leitura do disco usada pelos processos
* @param tbloc variável específica do algoritmo bCubing
* @param output_folder arquivo de saída, se o cubo for escrito em disco
* @param queries as consultas que serão executadas
* @param cube_algorithm o algoritmo de cubo de dados que será executado
* @param cube_table o nome de arquivo para o dataset
* @param on_demand booleano que determina se a computação será feita "on-demand"
*
*/
bool Handler::ParseInput(int argc, char *argv[], int my_rank,
                int num_procs, int &num_dims, int &num_meas, int &num_tuples,
                int &tuple_partition_size, int &dim_partition_size, int &reading_rate, int &tbloc, std::string &output_folder,
                std::vector<std::vector<int>> &queries,
                std::string &cube_algorithm, std::string &cube_table, bool &on_demand)
{
        try
        {
        		//Arquivo opcional de entrada com configurações do programas (as mesmas que também podem ser passadas por linha de comando)
                std::string user_config;

                //Variável que armazena o algoritmo escolhido pelo usuário
                std::string user_algorithm;

                //Variável que armazena o tipo de particionamento escolhido pelo usuário
                std::string user_partition;

                //Variável que armazena o nome do arquivo com tabela de entrada do algoritmo
                std::string user_table;

                //Caminho do sistema de arquivos onde o cubo deve ser escrito se necessário
                boost::filesystem::path user_output;

                //Lista de consultas do usuário ao cubo de dados
                std::vector<std::string> user_queries;

                //Lista de operações do usuário para executar em medidas cubo de dados
                std::vector<std::string> user_queries_ops;

                //Lista de operações conhecidas, aceitas como entradas válidas
                std::vector<std::string> known_ops { "soma", "Sm", "max", "Ma", "min",
                                "Mn", "média", "Me", "variância", "Va", "desvio padrão", "Dp",
                                "mediana", "Md", "moda", "Mo", "no-op", "No" };

                //Lista de opções do programa
                po::options_description generic("Opções gerais do programa");

                //Adicionando as opções à lista de opções (cada uma é um parâmetro disponível)
                generic.add_options()("dimensions,d", po::value<int>(),
                                "número absoluto de dimensões do cubo de dados")(
                                "measures,m", po::value<int>(),
                                "número absoluto de medidas do cubo de dados")(
                                "tuples,t", po::value<int>(),
                                "número absoluto de tuplas da tabela de entrada")(
                                "reading-rate,r", po::value<int>(),
                                "tarefas leem tuplas do disco em segmentos de tamanho n/r, onde n=t/processes")(
                                "tbloc",po::value<int>()->default_value(2),"tamanho de bloco (apenas algoritmo bcubing)")(
                                "output-folder,o",po::value<boost::filesystem::path>(&user_output),"um diretório em disco para escrever o cubo após computado")(
                                "query",
                                po::value<std::vector<std::string>>(
                                                &user_queries)->multitoken(),
                                "lista de consultas em formato string, entre aspas, separadas por espaço")(
                                "query-ops",
                                po::value<std::vector<std::string>>(
                                                                &user_queries_ops)->multitoken(),
                                "lista de operações (soma Sm, max Ma, min Mn, média Me, variância Va, desvio padrão Dp, mediana Md, moda Mo ou no-op No), uma por medida, em formato string, entre aspas, separadas por espaço")(
                                "algorithm,alg",
                                po::value<std::string>(&user_algorithm)->default_value(
                                                "bcubing"),
                                "tipo de cubo a ser criado (bcubing,fragcubing)")(
								"partition",
								po::value<std::string>(&user_partition)->default_value(
												"tuples"),
								"tipo de particionamento a ser usado (tuples,dimensions)")(
                                "table,f",
                                po::value<std::string>(&user_table),
                                "arquivo binário contendo uma tabela com dados para os quais deseja computar o cubo")(
                                "on-demand",
                                po::bool_switch(&on_demand)->default_value(
                                                false), "se true, cubo será computado com base nas consultas (mais rápido)");

                po::options_description config("Opção de arquivo de configuração");
                config.add_options()("config-file",
                                po::value<std::string>(&user_config),
                                "arquivo contendo uma ou mais das opções gerais, uma por linha");

                po::options_description other("Outras opções");
                other.add_options()("help", "mostra essa mensagem com as opções do programa")("version",
                                "imprime a versão do programa");

                //Por linha de comando pode usar qualquer uma das opções
                po::options_description command_line;
                command_line.add(generic).add(config).add(other);

                //No arquivo de configuração deve passar apenas alguma das opções gerais
                po::options_description command_file;
                command_file.add(generic);

                //Um mapa de variáveis, guarda as variáveis passadas como opções do programa
                po::variables_map vm;

                //Processa as entradas por linha de comando
                po::store(
                                po::command_line_parser(argc, argv).options(
                                                command_line).run(), vm);

                //Aqui de fato atualiza o mapa de variáveis
                po::notify(vm);

                //Se usuário digitou --help ou não passou nenhuma opção, mostre um resumo das opções
                if (vm.count("help") || (argc <= 1))
                {
                        std::cout
                                        << "datacube-mpi-framework - um facilitador para cubos de dados em cluster"
                                        << std::endl << std::endl;
                        std::cout
                                        << "Utilização: mpirun datacube-mpi-framework [options] [table-file]"
                                        << std::endl
										<< "Exemplo: mpirun -np 2 datacube-mpi-framework -d 4 -m 1 -t 10 -r 1 --query \"? * 1 *\" -f dataset_10t_4d.bin -o cubo/ --algorithm bcubing --tbloc 4"
										<< std::endl;

                        //Resuno das opções por linha de comando
                        std::cout << command_line << std::endl;
                        return false;
                }

                //Se digitou --version mostre a versão atual do programa
                if (vm.count("version"))
                {
                        std::cout << "datacube-mpi-framework, v1.0"
                                        << std::endl;
                        return false;
                }

                //Se digitou --config-file, tente abrir e ler o arquivo de configuração
                if (vm.count("config-file"))
                {
                		//Abertura do arquivo de configuração
                        std::ifstream ifconf(user_config.c_str());

                        if (!ifconf)
                        {
                                std::cout << "ERRO: Não é possível abrir o arquivo "
                                                << user_config << std::endl;
                                return false;
                        }
                        else
                        {
                        		//Armazena os parâmetros do arquivo (um por linha)
                                po::store(
                                                parse_config_file(ifconf,
                                                                command_file),
                                                vm);

                                //Atualiza o mapa de variáveis
                                notify(vm);
                        }
                }

                //Se digitou --output-folder, verifique se o diretório existe ou pode ser criado
                if (vm.count("output-folder")){
                	//Se o diretório não existe
                	if(!boost::filesystem::exists(user_output)){
                		//Apenas o primeiro processo se encarrega de tentar criar o diretório
                		if(my_rank == 0){
                    		//Tente criar, mas se não conseguir informe erro ao usuário
                    		if(!boost::filesystem::create_directory(user_output)){
        						std::cout << "ERRO: O diretório "
        										<< user_output
        										<< " não existe e não pode ser criado." << std::endl;
        						return false;
                    		}
                		}
                	}
                }
                else{
                	//O algoritmo bcubing exige a existência de um diretório de saída
                	if(user_algorithm == "bcubing"){
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "output-folder");
                        return false;
                	}
                }

                //É obrigatório passar uma o caminho de uma tabela com dados de entrada
                if (!vm.count("table"))
                {
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "table");
                        return false;
                }
                else
                {
                		//Tenta abrir a tabela indicada pelo usuário
                        std::ifstream iftable(user_table.c_str());

                        if (!iftable)
                        {
                                std::cout << "ERRO: Não é possível abrir a tabela "
                                                << user_table << std::endl;
                                return false;
                        }
                }

                //É obrigatório passar a quantidade de dimensões da tabela
                if (!vm.count("dimensions"))
                {
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "dimensions");
                        return false;
                }

                //É obrigatório passar a quantidade de medidas da tabela
                if (!vm.count("measures"))
                {
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "measures");
                        return false;
                }

                //É obrigatório passar a quantidade de tuplas da tabela
                if (!vm.count("tuples"))
                {
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "tuples");
                        return false;
                }

                //É obrigatório informar a taxa de leitura da tabela do disco
                if (!vm.count("reading-rate"))
                {
                        throw po::validation_error(
                                        po::validation_error::at_least_one_value_required,
                                        "reading-rate");
                        return false;
                }

                //É obrigatório informar o algoritmo a ser utilizado
                if (vm.count("algorithm")
                                && (!(user_algorithm == "bcubing"
                                                || user_algorithm
                                                                == "fragcubing")))
                {
                        throw po::validation_error(
                                        po::validation_error::invalid_option_value,
                                        "algorithm");
                        return false;
                }

                //Salva o número de dimensões com base no valor passado pelo usuário
                num_dims = std::abs(vm["dimensions"].as<int>());

                //Salva o número de medias com base no valor passado pelo usuário
                num_meas = std::abs(vm["measures"].as<int>());

                //Salva o número de tuplas com base no valor passado pelo usuário
                num_tuples = std::abs(vm["tuples"].as<int>());

                //Taxa de leitura de dados do disco (para ler a tabela de entrada)
                reading_rate = std::abs(vm["reading-rate"].as<int>());

                //Tamanho do bloco, aplicável somente ao algoritmo bCubing
                tbloc = std::abs(vm["tbloc"].as<int>());

                //Salva o algoritmo de cubos escolhido
                cube_algorithm = user_algorithm;

                //Salva a tabela de dados a ser usada para computar o cubo
                cube_table = user_table;

                //Salva a informação do diretório de saída para o cubo de dados
                output_folder = boost::filesystem::weakly_canonical(user_output).string();

                //Verifica se o usuário informou alguma query para ser executada e já tenta converter para inteiros
                if (vm.count("query")
                                && !(QueryConverter(user_queries, queries,
                                                num_dims)))
                {
                        std::cout
                                        << "ERRO: por favor garanta que todas as consultas são válidas antes de tentar novamente."
                                        << std::endl;
                        return false;
                }

                //Não é obrigatório, mas o usuário pode fornecer uma lista de operações para cada medida (frequência é sempre calculado)
                if (vm.count("query-ops")) {

                        //Se passar mais listas de operações do que consultas, significa que pelo menos uma lista de operações é desnecessária
                        if (user_queries_ops.size() > user_queries.size()) {
                                std::cout
                                                << "ERRO: você passou mais listas de operações do que consultas a serem executadas."
                                                << std::endl;
                                return false;
                        }

                        //Agora vamos verificar cada lista de operações fornecida
                        for (auto &operations : user_queries_ops) {

                                //Será usado para contar quantos operadores foram fornecidos
                                int num_args = 0;

                                //Armazena um dos operadores fornecidos
                                std::string arg;

                                //A lista de operações, em formato de string
                                std::stringstream ops(operations);

                                //Passa por cada argumento da lista de operações, separado por espaços
                                while (ops >> arg) {
                                        //Apenas alguns operadores são aceitos, se passar algum errado deve encerrar a execução
                                        if (std::find(std::begin(known_ops), std::end(known_ops), arg) == std::end(known_ops)) {
                                                std::cout
                                                                << "ERRO: o operador " << arg << " na lista de operações \"" << operations << "\" não é válido."
                                                                << std::endl;
                                                return false;
                                        }
                                        num_args++;
                                }

                                //Se chegou aqui é porque todos os operadores são válidos, mas talvez tenha passado operadores demais/menos
                                if(num_args != num_meas){
                                        std::cout
                                                        << "ERRO: o número de operadores na lista de operações \"" << operations << "\" é diferente do número de medidas existentes."
                                                        << std::endl;
                                        return false;
                                }
                        }
                }

                //Tenta obter o tamanho da partição dos dados pelo número de processos em execução
                if (!(GetPartitionSize(num_tuples, num_dims, num_procs, my_rank,
                		tuple_partition_size, dim_partition_size, reading_rate, user_partition)))
                {
                        std::cout
                                        << "ERRO: por favor garanta que os valores de tuplas, tarefas e a taxa de leitura sejam corretos."
                                        << std::endl;
                        return false;

                }

                //Mensagem final para o usuário apenas leitura dos parâmetros
                if(tuple_partition_size > 0){
                    std::cout << "I am rank " << my_rank << " and will take a partition of " << tuple_partition_size << " tuples with " << num_dims << " dimensions..." << std::endl;
                } else if (dim_partition_size > 0){
                    std::cout << "I am rank " << my_rank << " and will take a partition of " << num_tuples << " tuples with " << dim_partition_size << " dimensions..." << std::endl;
                }


        }
        catch (std::exception& e)
        {
                std::cout << "ERRO: " << e.what() << std::endl;
                return false;
        }
        return true;
}
