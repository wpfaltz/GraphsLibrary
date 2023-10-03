#include "../library/GraphLibrary.h"
#include <random>
#include <chrono>
#include <thread>

using namespace std;

int graphs[1] = {4};

int main() {
    for (int j = 0; j < sizeof(graphs) / sizeof(graphs[0]); j++) {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        std::cout << "Temporizador acabado, iniciando leitura" << std::endl;
        AdjacencyMatrixGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[j]);
        std::string pathToGraph = "./graphs/grafo_" + graphNumberToString + ".txt";
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
        std::cout << "Leitura finalizada, aguardando 15 segundos" << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(15));
    }
    return 0;
}
