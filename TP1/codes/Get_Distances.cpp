#include "../../library/GraphLibrary.h"
#include <random>
#include <chrono>

using namespace std;

int graphs[6] = {1,2,3,4,5,6};

int main() {
        for (int j = 0; j < sizeof(graphs) / sizeof(graphs[0]); j++) {
                AdjacencyListGraph graph(0); // Crie o grafo com 0 vértices
                std::string graphNumberToString = std::to_string(graphs[j]);
                std::string pathToGraph = "./graphs/grafo_" + graphNumberToString + ".txt";
                graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo
                std::cout << "Grafo " << graphNumberToString << ":" << std::endl;
                int distancia1 = graph.getDistance(10, 20);
                std::cout << "  A distância entre os vértices 10 e 20 é " << distancia1 << std::endl;
                int distancia2 = graph.getDistance(10, 30);
                std::cout << "  A distância entre os vértices 10 e 30 é " << distancia2 << std::endl;
                int distancia3 = graph.getDistance(20, 30);
                std::cout << "  A distância entre os vértices 20 e 30 é " << distancia3 << std::endl;
        }
}
