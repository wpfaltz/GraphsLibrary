#include "../../library/GraphLibrary.h"
#include <random>
#include <chrono>

using namespace std;

int graphs[1] = {1};

int main() {
    for (int j = 0; j < sizeof(graphs) / sizeof(graphs[0]); j++) {
        AdjacencyListGraph graph(0); // Crie o grafo com 0 vértices
        std::string graphNumberToString = std::to_string(graphs[j]);
        std::string pathToGraph = "./graphs/grafo_" + graphNumberToString + ".txt";
        graph.readGraphFromFile(pathToGraph); // Leia o grafo do arquivo

        std::string pathToGraphInfo = "./outputs/Info/Grafo_" + graphNumberToString + "_Info.txt";
        graph.writeGraphInfoToFile(pathToGraphInfo);

        int numVertices = graph.getNumVertices();

        int startVertices[1] = {};
        for (int i = 0; i < sizeof(startVertices) / sizeof(startVertices[0]); i++) {
            std::string numberToString = std::to_string(startVertices[i]);
            std::string pathToBFSOutput = "./outputs/Search_Results/Graph_" + graphNumberToString + "/BFS/Vertex_" + numberToString + ".txt";
            graph.BFSWithOutput(startVertices[i], pathToBFSOutput);
            std::string pathToDFSOutput = "./outputs/Search_Results/Graph_" + graphNumberToString + "/DFS/Vertex_" + numberToString + ".txt";
            graph.DFSWithOutput(startVertices[i], pathToDFSOutput);
        }
        
        
    }
    
    return 0;
}
