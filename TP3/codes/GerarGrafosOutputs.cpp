#include "../library/WeightedGraphs.h"
#include <random>
#include <chrono>
#include <thread>

#include <bits/stdc++.h> 
using namespace std;

int main() {
    for (int i=1; i<7; i++) {
        WeightedGraph graph(0, true); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(i);
        std::string pathToGraph = "./TP3/graphs/grafo_rf_" + graphNumberToString + ".txt";

        std::cout << "- Grafo " + graphNumberToString + ":" << std::endl;
        
        std::cout << "Inicializando leitura do arquivo correspondente ao grafo" << std::endl;
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        std::cout << "Leitura concluída com sucesso, o grafo está representado em memória!" << std::endl;

        std::string pathToOutput = "./TP3/outputs/grafo_rf_" + graphNumberToString + ".txt";
        graph.Output_Ford_Fulkerson_BFS(pathToOutput);

        std::cout << "Inicializando execução que fará a escrita no arquivo" << std::endl;
        auto result = graph.Ford_Fulkerson_BFS(1, 2, true);
        std::cout << "Grafo residual escrito para o arquivo com sucesso." << std::endl;
        std::cout << "" << std::endl;
    }
}
