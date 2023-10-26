#include "../../library/WeightedGraphs.h"
#include <random>
#include <chrono>
#include <iomanip>

#include <bits/stdc++.h> 
using namespace std;

int graphs[1] = {1};

int main() {
    for (int i = 0; i < sizeof(graphs) / sizeof(graphs[0]); i++) {
        WeightedGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[i]);
        std::string pathToGraph = "./graphs/TP2/rede_colaboracao.txt";
        std::cout << ". rede_colaboracao:" << std::endl;
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        std::cout << " Grafo lido com sucesso! " << std::endl;

        int startVertex = 2722;
        int finalVertices[5] =  {11365, 471365, 5709, 11386, 343930};

        auto start = std::chrono::high_resolution_clock::now();
        graph.DijkstraWithHeap(startVertex);
        auto end = std::chrono::high_resolution_clock::now();
        double executionTimeHeap = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
        std::cout << "  Tempo de execução com Heap: " << executionTimeHeap << " segundos" << std::endl;

        for (int j = 0; j < sizeof(finalVertices) / sizeof(finalVertices[0]); j++) {
            std::cout << "      A menor distância do vértice " << startVertex << " ao vértice " << finalVertices[j] << " é " << graph.getMinDistance(finalVertices[j]) << std::endl;

            std::list<int> minimalPath = graph.minimalPathToVertex(finalVertices[j]);

            std::cout << "      O caminho mínimo do vértice " << startVertex << " ao vértice " << finalVertices[j] << " é dado por: ";
            for (auto it = minimalPath.begin(); it != minimalPath.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
    }
    return 0;
}
