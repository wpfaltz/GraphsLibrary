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

        std::cout << "Iniciando as 10 execuções de teste" << std::endl;
        int rodadas = 10;
        float total_execution_time = 0.0;
        for (int i = 0; i < rodadas+1; i++) {
            if (i == 0) {
                auto start = std::chrono::high_resolution_clock::now();
                auto result = graph.Ford_Fulkerson_BFS(1, 2, false);
                auto end = std::chrono::high_resolution_clock::now();
                double executionTimeFord_Fulkerson_BFS = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
                total_execution_time += executionTimeFord_Fulkerson_BFS;
                int maxFlow = result.first;
                std::cout << "Fluxo máximo: " << maxFlow << std::endl;
            }

            else if (i == 10) {
                std::cout << "10 execuções para aferição do tempo médio de execução finalizadas" << std::endl;
                float mean_execution_time = total_execution_time / rodadas;
                std::cout << "Tempo médio de execução do Ford-Fulkerson: " << mean_execution_time << " segundos" << std::endl;
                std::cout << "Inicializando execução que fará a escrita no arquivo" << std::endl;
                auto result = graph.Ford_Fulkerson_BFS(1, 2, true);
                std::cout << "Grafo residual escrito para o arquivo com sucesso." << std::endl;
            }

            else {
                auto start = std::chrono::high_resolution_clock::now();
                auto result = graph.Ford_Fulkerson_BFS(1, 2, false);
                auto end = std::chrono::high_resolution_clock::now();
                double executionTimeFord_Fulkerson_BFS = std::chrono::duration<double, std::nano>(end - start).count() / 1e9;
                total_execution_time += executionTimeFord_Fulkerson_BFS;
            }
        }
        std::cout << "" << std::endl;
    }
}
