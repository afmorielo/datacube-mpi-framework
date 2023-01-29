#include "DataCube.h"
#include <algorithm>
#include <cassert>

// Aqui posso declarar uma única vez esssa variáveis
// globais e qualquer um que inclua o header pode utilizá-las
MPI_Datatype MPI_TUPLE_DIMS, MPI_TUPLE_MEAS;

DataCube::~DataCube()
{

}

std::vector<int> DataCube::IntersectTwoVectors(std::vector<int> const& vector_a, std::vector<int> const& vector_b) {

	//Primeira posição da primeira lista de inteiros
    auto a = vector_a.begin();

    //Posição final da primeira lista de inteiros
    auto a_end = vector_a.end();

	//Primeira posição da segunda lista de inteiros
    auto b = vector_b.begin();

    //Posição final da segunda lista de inteiros
    auto b_end = vector_b.end();

    //Verifica se a primeira lista está ordenada (necessário para o algoritmo de interseção)
    assert(std::is_sorted(a, a_end));

    //Verifica se a segunda lista está ordenada (necessário para o algoritmo de interseção)
    assert(std::is_sorted(b, b_end));

	//Irá guardar a interseção final das múltiplas listas de inteiros
    std::vector<int> intersection_vector;

    //Percorre as duas listas de inteiros simultaneamente
    while (a != a_end && b != b_end) {

    	//Se encontrou valores iguais, adicione à interseção e incremente os iteradores
        if (*a == *b) {
        	intersection_vector.push_back(*a);
            ++a;
            ++b;

            //Força a pular para a próxima etapa do laço
            continue;
        }

        //Não são iguais, se o primeiro for menor, incremente seu iterador
        if (*a < *b) {
            ++a;

            //Força a pular para a próxima etapa do laço
            continue;
        }

        //Teste de lógica: se não são iguais e o primeiro não é menor então o segundo deve ser menor
        assert(*a > *b);

        //Incremente o iterador do segundo e vá para a próxima etapa do laço
        ++b;
    }

	//Retorna a interseção final das listas de inteiros
    return intersection_vector;
}

std::vector<int> DataCube::IntersectMultipleVectors(std::vector<std::vector<int>> &sorted_vectors){

	//Irá guardar a interseção final das múltiplas listas de inteiros
	//O algortimo inicialmente define a interseção como a primeira lista
	std::vector<int> intersection_vector = sorted_vectors[0];

	//Armazena vetores intermediários, interseções geradas em cada iteração do algoritmo
	std::vector<int> intermediate;

	//Começa na posição 1, fazendo a interseção dois a das listas de números
	for (std::vector<T>::size_type vector_index = 1; vector_index < sorted_vectors.size(); ++vector_index) {

		//A interseção intermediária é resultado da interseção atual com a lista seguinte de inteiros
		intermediate = IntersectTwoVectors(intersection_vector, sorted_vectors[vector_index]);

		//Troca o valor da interseção atual com a interseção intermediária
		std::swap(intersection_vector, intermediate);

		//Limpa o espaço em memória da interseção intermediária
		intermediate.clear();
	}

	//Retorna a interseção final das listas de inteiros
	return intersection_vector;
}
