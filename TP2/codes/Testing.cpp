#include "../library/WeightedGraphs.h"
#include <random>
#include <chrono>
#include <iomanip>

#include <bits/stdc++.h> 
using namespace std;

int graphs[1] = {4};

int main() {
    for (int i = 0; i < sizeof(graphs) / sizeof(graphs[0]); i++) {
        WeightedGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[i]);
        std::string pathToGraph = "./TP2/graphs/grafo_W_" + graphNumberToString + ".txt";
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        std::cout << "Grafo " << graphs[i] << " lido com sucesso! " << std::endl;
        int numVertices = graph.getNumVertices();

        // Defina o intervalo de números possíveis para o vértice inicial aleatório
        int minNumber = 1;
        int maxNumber = numVertices;

        // Crie um objeto de geração de números aleatórios
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distribution(minNumber, maxNumber);

        // Salve os números aleatórios em um arquivo
        std::string pathToChosenVertices = "./TP2/outputs/grafo_W_" + graphNumberToString + "/Chosen_Vertices.txt";
        std::ofstream outputVerticesFile(pathToChosenVertices);
        for (int i = 0; i < 100; ++i) {
            int randomNumber = distribution(gen);
            outputVerticesFile << randomNumber << std::endl;
        }
        outputVerticesFile.close();
        std::cout << "Vértices aleatórios escolhidos escritos no arquivo com sucesso!" << std::endl;

        std::string pathToElapsedTimesHeap = "./TP2/outputs/grafo_W_" + graphNumberToString + "/Heap_Elapsed_Times.txt";
        std::ofstream elapsedTimesHeapFile(pathToElapsedTimesHeap);

        std::string pathToElapsedTimesVector = "./TP2/outputs/grafo_W_" + graphNumberToString + "/Vector_Elapsed_Times.txt";
        std::ofstream elapsedTimesVectorFile(pathToElapsedTimesVector);

        std::ifstream inputVerticesFile(pathToChosenVertices);
        if (inputVerticesFile.is_open()) {
            int readVertex;
            while (inputVerticesFile >> readVertex) {
                auto start = std::chrono::high_resolution_clock::now();
                graph.DijkstraWithHeap(readVertex);
                auto end = std::chrono::high_resolution_clock::now();
                double executionTimeHeap = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
                elapsedTimesHeapFile << executionTimeHeap << std::endl;

                start = std::chrono::high_resolution_clock::now();
                graph.DijkstraWithVector(readVertex);
                end = std::chrono::high_resolution_clock::now();
                double executionTimeVector = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
                elapsedTimesVectorFile << executionTimeVector << std::endl;
            }
        }
        elapsedTimesHeapFile.close();
        elapsedTimesVectorFile.close();
        inputVerticesFile.close();
    }
    return 0;
}
