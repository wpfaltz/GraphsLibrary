#include "../../library/GraphLibrary.h"
#include <random>
#include <chrono>

using namespace std;

int graphs[6] = {1,2,3,4,5,6};

int main() {
        for (int i = 0; i < sizeof(graphs) / sizeof(graphs[0]); i++) {
                AdjacencyListGraph graph(0); // Crie o grafo com 0 vértices
                std::string graphNumberToString = std::to_string(graphs[i]);
                std::string pathToGraph = "./graphs/grafo_" + graphNumberToString + ".txt";
                graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
                std::cout << "Grafo " << graphNumberToString << ":" << std::endl;
                int diameter = graph.getDiameter();
                std::cout << "  O diâmetro do grafo é " << diameter << std::endl;
        
        }
}
