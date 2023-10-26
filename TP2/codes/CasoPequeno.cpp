#include "../../library/WeightedGraphs.h"
#include <random>
#include <chrono>
#include <thread>

#include <bits/stdc++.h> 
using namespace std;

int main() {
    WeightedGraph graph(0); // Crie o grafo com 0 vértices
    std::string pathToGraph = "./graphs/TP2/rede_colaboracao.txt";
    graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo

    graph.DijkstraWithHeap(2722);
    std::cout << "  Dijkstra com vetor: " << std::endl;
    std::cout << "      A distância mínima de Edsger W. Dijkstra à Alan Turing é " << graph.getMinDistance(11365) << std::endl;
    std::cout << "      A distância mínima de Edsger W. Dijkstra à J. B. Kruskal é " << graph.getMinDistance(471365) << std::endl;
    std::cout << "      A distância mínima de Edsger W. Dijkstra à Jon M. Kleinberg é " << graph.getMinDistance(5709) << std::endl;
    std::cout << "      A distância mínima de Edsger W. Dijkstra à Éva Tardos é " << graph.getMinDistance(11386) << std::endl;
    std::cout << "      A distância mínima de Edsger W. Dijkstra à Daniel R. Figueiredo é " << graph.getMinDistance(343930) << std::endl;
        
    return 0;
}
