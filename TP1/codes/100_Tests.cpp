#include "../../library/GraphLibrary.h"
#include <random>
#include <chrono>

using namespace std;

int graphs[1] = {6};

int main() {
    for (int j = 0; j < sizeof(graphs) / sizeof(graphs[0]); j++) {
        AdjacencyListGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[j]);
        std::string pathToGraph = "./graphs/grafo_" + graphNumberToString + ".txt";
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        int numVertices = graph.getNumVertices();

        // Defina o intervalo de números possíveis
        int minNumber = 1;
        int maxNumber = numVertices;

        // Crie um objeto de geração de números aleatórios
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distribution(minNumber, maxNumber);

        // Salve os números aleatórios em um arquivo
        std::string pathToChosenNumbers = "./outputs/Search_Results/Graph_" + graphNumberToString + "/Chosen_Numbers.txt";
        std::ofstream outputNumbersFile(pathToChosenNumbers);
        for (int i = 0; i < 100; ++i) {
            int randomNumber = distribution(gen);
            outputNumbersFile << randomNumber << std::endl;
        }
        outputNumbersFile.close();

        std::string pathToElapsedTimesBFS = "./outputs/Search_Results/Graph_" + graphNumberToString + "/BFS_Elapsed_Times.txt";
        std::ofstream elapsedTimesBFSFile(pathToElapsedTimesBFS);

        std::string pathToElapsedTimesDFS = "./outputs/Search_Results/Graph_" + graphNumberToString + "/DFS_Elapsed_Times.txt";
        std::ofstream elapsedTimesDFSFile(pathToElapsedTimesDFS);
        
        std::ifstream inputNumbersFile(pathToChosenNumbers);
        if (inputNumbersFile.is_open()) {
            int readNumber;
            while (inputNumbersFile >> readNumber) {
                auto startTimeBFS = std::chrono::high_resolution_clock::now();
                graph.BFSNoOutput(readNumber);
                auto endTimeBFS = std::chrono::high_resolution_clock::now();
                auto durationBFS = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeBFS - startTimeBFS);
                long long executionTimeBFS = durationBFS.count();
                elapsedTimesBFSFile << executionTimeBFS << std::endl;

                auto startTimeDFS = std::chrono::high_resolution_clock::now();
                graph.DFSNoOutput(readNumber);
                auto endTimeDFS = std::chrono::high_resolution_clock::now();
                auto durationDFS = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDFS - startTimeDFS);
                long long executionTimeDFS = durationDFS.count();
                elapsedTimesDFSFile << executionTimeDFS << std::endl;

            }
            inputNumbersFile.close();
        }
    }

    return 0;
}
