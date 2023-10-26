#include "../../library/WeightedGraphs.h"
#include <random>
#include <chrono>
#include <iomanip>

#include <bits/stdc++.h> 
using namespace std;

int graphs[5] = {1,2,3,4,5};

int main() {
    for (int i = 0; i < sizeof(graphs) / sizeof(graphs[0]); i++) {
        WeightedGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[i]);
        std::string pathToGraph = "./graphs/TP2/grafo_W_" + graphNumberToString + ".txt";
        std::cout << ". grafo_W_" + graphNumberToString + ":" << std::endl;
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        std::cout << " Grafo " << graphs[i] << " lido com sucesso! " << std::endl;

        int startVertex = 10;
        int finalVertices[5] =  {20, 30, 40, 50, 60};

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
        

        start = std::chrono::high_resolution_clock::now();
        graph.DijkstraWithVector(startVertex);
        end = std::chrono::high_resolution_clock::now();
        double executionTimeVector = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
        std::cout << "Tempo de execução com vetor: " << executionTimeVector << " segundos" << std::endl;

        for (int j = 0; j < sizeof(finalVertices) / sizeof(finalVertices[0]); j++) {
            std::cout << "A menor distância do vértice " << startVertex << " ao vértice " << finalVertices[j] << " é " << graph.getMinDistance(finalVertices[j]) << std::endl;

            std::list<int> minimalPath = graph.minimalPathToVertex(finalVertices[j]);

            std::cout << "O caminho mínimo do vértice " << startVertex << " ao vértice " << finalVertices[j] << " é dado por: ";
            for (auto it = minimalPath.begin(); it != minimalPath.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}
