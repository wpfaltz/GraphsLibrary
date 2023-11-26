#ifndef WEIGHTEDGRAPHS_H
#define WEIGHTEDGRAPHS_H

#include <bits/stdc++.h> 
using namespace std;
#include <iostream>

class WeightedGraph {
public:
    struct Edge {
        int capacity;
        int flow;
    };

    // Construtor que define o número de vértices
    WeightedGraph(int numVertices, bool directedGraph);

    // Método para construir o grafo a partir do arquivo lido
    void readGraphFromFile(const std::string& filename);

    // Métodos para adicionar e remover arestas
    void addEdge(int vertexA, int vertexB, int capacityAB);
    void removeEdge(int vertexA, int vertexB);

    // Demais métodos
    int getNumVertices();
    int getNumEdges();
    void Output_Ford_Fulkerson_BFS(const std::string filename);
    std::pair<int, std::vector<std::map<int, WeightedGraph::Edge>>> Ford_Fulkerson_BFS(int source, int sink, bool saveToDisk);

    // void Ford_Fulkerson_Scaled(int source, int sink);

private:
    bool hasNegativeDistance = false;
    int numVertices;
    int numEdges;
    int numComponents;
    bool directedGraph;
    std::vector<std::map<int, Edge>> edgeMaps;
    std::vector<std::map<int, int>> residualGraph;
    std::vector<int> dist; // Vetor de distâncias
    std::vector<int> predecessors; // vetor de predecessores
    std::string output_file_path_ford_fulkerson_bfs;

    std::list<int> findAugmentingPathBFS(std::vector<std::map<int, int>>& residualGraph, int source, int sink);
    void constructInitialResidualGraph(const std::vector<std::map<int, Edge>>& edgeMaps);
    std::pair<std::pair<int, int>, int> findBottleneck(const std::list<int>& path, const std::vector<std::map<int, int>>& residualGraph);
    void updateResidualGraph(std::vector<std::map<int, int>>& residualGraph, const std::list<int>& path, int bottleneckFlow, std::pair<int, int> bottleneckEdge);
    void updateOriginalGraph(const std::list<int>&path, int bottleneckFlow, std::pair<int, int> bottleneckEdge);
    void saveGraphToFile(const std::string& filename, int maxFlow);
};

WeightedGraph::WeightedGraph(int numVertices, bool directedGraph) {
    this->numVertices = numVertices;
    this->numEdges = 0;
    this->directedGraph = directedGraph;
    this->dist.resize(numVertices, std::numeric_limits<int>::infinity()); // Inicializa o vetor de distâncias
    this->predecessors.resize(numVertices, -1);
    this->edgeMaps.resize(numVertices);
}

void WeightedGraph::addEdge(int vertexA, int vertexB, int capacityAB) {
    Edge edgeAB = {capacityAB, 0};
    edgeMaps[vertexA - 1][vertexB] = edgeAB;
    numEdges++;
}

void WeightedGraph::removeEdge(int vertexA, int vertexB) {
    auto itA = edgeMaps[vertexA - 1].find(vertexB);
    if (itA != edgeMaps[vertexA - 1].end()) {
        edgeMaps[vertexA - 1].erase(itA);
    }

    // Se o grafo não for direcionado, remover a aresta reversa também
    if (!directedGraph) {
        auto itB = edgeMaps[vertexB - 1].find(vertexA);
        if (itB != edgeMaps[vertexB - 1].end()) {
            edgeMaps[vertexB - 1].erase(itB);
        }
    }
}

void WeightedGraph::readGraphFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    inputFile >> numVertices;
    
    edgeMaps.resize(numVertices);
    dist.assign(numVertices, std::numeric_limits<int>::infinity()); // Inicializa o vetor de distâncias para cada vértice

    int vertexA, vertexB ;
    int capacityAB;
    int edgeCount = 0;
    while (inputFile >> vertexA >> vertexB >> capacityAB) {
        if (capacityAB < 0) {
            hasNegativeDistance = true;
        }

        addEdge(vertexA, vertexB, capacityAB);

        if (directedGraph == false) {
            addEdge(vertexB, vertexA, capacityAB);
        }
    }

    inputFile.close();
}

int WeightedGraph::getNumEdges(){
    if (directedGraph == false) {
        return numEdges/2;
    }
    else {
        return numEdges;
    }
}

int WeightedGraph::getNumVertices(){
    return numVertices;
}

void WeightedGraph::Output_Ford_Fulkerson_BFS(const std::string filename){
    output_file_path_ford_fulkerson_bfs = filename;
}

void WeightedGraph::constructInitialResidualGraph(const std::vector<std::map<int, Edge>>& edgeMaps) {
    // Certificando que o residualGraph está vazio antes de começar
    residualGraph.clear();

    residualGraph.resize(edgeMaps.size());

    // Para cada vértice, construa o mapa de destinos associados a arestas no residualGraph
    for (size_t i = 0; i < edgeMaps.size(); ++i) {
        // Mapa correspondente a esta posição no residualGraph
        std::map<int, int>& residualEdges = residualGraph[i];

        for (const auto& entry : edgeMaps[i]) {
            int destination = entry.first;
            const Edge& edge = entry.second;

            // Adicionar a aresta no residualGraph
            residualEdges[destination] = edge.capacity;
        }
    }
}

std::pair<int, std::vector<std::map<int, WeightedGraph::Edge>>> WeightedGraph::Ford_Fulkerson_BFS(int source, int sink, bool saveToDisk) {

    // Inicializar o grafo residual
    constructInitialResidualGraph(edgeMaps);

    // Inicializar o fluxo máximo
    int maxFlow = 0;

    // Continuar até não haver mais caminhos aumentantes
    while (true) {
        // Encontrar um caminho aumentante usando BFS
        std::list<int> path = findAugmentingPathBFS(residualGraph, source, sink);

        // Se não houver mais caminhos aumentantes, sair do loop
        if (path.empty()) {
            break;
        }

        // Encontrar a capacidade mínima ao longo do caminho aumentante
        int bottleneckFlow;
        std::pair<int, int> bottleneckEdge;
        std::tie(bottleneckEdge, bottleneckFlow) = findBottleneck(path, residualGraph);


        // Atualizar a rede residual com o fluxo encontrado
        updateResidualGraph(residualGraph, path, bottleneckFlow, bottleneckEdge);

        // Atualizar o grafo original com o fluxo
        updateOriginalGraph(path, bottleneckFlow, bottleneckEdge);

        // Adicionar o fluxo ao fluxo máximo
        maxFlow += bottleneckFlow;
    }

    
    if (saveToDisk == true) {
        saveGraphToFile(output_file_path_ford_fulkerson_bfs, maxFlow);
    }

    // O fluxo máximo está armazenado em maxFlow
    // A rede residual está armazenada em residualGraph
    return {maxFlow, edgeMaps};
}

std::list<int> WeightedGraph::findAugmentingPathBFS(std::vector<std::map<int, int>>& residualGraph, int source, int sink) {
    std::queue<int> queue;
    std::vector<int> parent(numVertices, -1);
    std::vector<bool> visited(numVertices, false);
    std::list<int> path;

    queue.push(source);
    visited[source - 1] = true;

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        for (const auto& entry : residualGraph[currentVertex - 1]) {
            int destination = entry.first;

            // Verificar se a aresta não está saturada e se o destino não foi visitado
            if (residualGraph[currentVertex - 1][destination] > 0 && !visited[destination - 1]) {
                // Armazenar o pai do destino
                parent[destination - 1] = currentVertex;

                // Marcar o destino como visitado
                visited[destination - 1] = true;

                // Se o destino for o sumidouro, reconstruir o caminho e encerrar a busca
                if (destination == sink) {
                    int current = destination;
                    while (current != source) {
                        path.push_front(current);  // Inserir no início da lista
                        current = parent[current - 1];
                    }
                    path.push_front(source);
                    return path;
                }

                // Adicionar o destino à fila para explorar seus vizinhos
                queue.push(destination);
            }
        }
    }

    // Se nenhum caminho aumentante for encontrado
    return path;
}

std::pair<std::pair<int, int>, int> WeightedGraph::findBottleneck(const std::list<int>& path, const std::vector<std::map<int, int>>& residualGraph) {
    int bottleneck = std::numeric_limits<int>::max();
    std::pair<int, int> bottleneckEdge = {-1, -1};

    auto it = path.begin();
    while (it != std::prev(path.end())) {
        int u = *it;
        int v = *std::next(it);

        // Encontrar a capacidade residual na aresta (u, v) do caminho
        auto itResidual = residualGraph[u - 1].find(v);
        if (itResidual != residualGraph[u - 1].end()) {
            int allocatedFlow = itResidual->second;

            // Atualizar o gargalo se a capacidade residual for menor
            if (allocatedFlow < bottleneck) {
                bottleneck = allocatedFlow;
                bottleneckEdge = {u, v};
            }
        }

        ++it;
    }

    return std::make_pair(bottleneckEdge, bottleneck);
}

void WeightedGraph::updateResidualGraph(std::vector<std::map<int, int>>& residualGraph, const std::list<int>& path, int bottleneckFlow, std::pair<int, int> bottleneckEdge) {

    for (auto it = path.begin(); it != path.end(); ++it) {
        int currentVertex = *it;
        auto nextIt = std::next(it);

        if (nextIt == path.end()) {
            break;  // Último vértice no caminho, não há aresta seguinte
        }

        int nextVertex = *nextIt;

        // Atualizar o fluxo na aresta (currentVertex, nextVertex)
        residualGraph[currentVertex - 1][nextVertex] -= bottleneckFlow;

        // Inicializar a aresta reversa se não existir
        if (residualGraph[nextVertex - 1].find(currentVertex) == residualGraph[nextVertex - 1].end()) {
            residualGraph[nextVertex - 1][currentVertex] = 0;  // Inicializando com capacidade zero
        }      

        // Atualizar o fluxo na aresta reversa (nextVertex, currentVertex)
        residualGraph[nextVertex - 1][currentVertex] += bottleneckFlow;

        // Remover aresta do grafo residual se a capacidade atingir zero
        if (residualGraph[currentVertex - 1][nextVertex] == 0) {
            residualGraph[currentVertex - 1].erase(nextVertex);
        }
    }
}

void WeightedGraph::updateOriginalGraph(const std::list<int>& path, int bottleneckFlow, std::pair<int, int> bottleneckEdge) {
    for (auto it = path.begin(); it != path.end(); ++it) {
        int currentVertex = *it;
        auto nextIt = std::next(it);

        if (nextIt == path.end()) {
            break;  // Último vértice no caminho, não há aresta seguinte
        }

        int nextVertex = *nextIt;

        if (edgeMaps[currentVertex - 1].find(nextVertex) == edgeMaps[currentVertex - 1].end()) {
            // Se nao encontrar, eh pq estamos procurando uma aresta reversa
            // Então parte do fluxo vai ser realocado, então diminui o fluxo que passa pela aresta original (nextVertex, currentVertex)
            edgeMaps[nextVertex - 1][currentVertex].flow -= bottleneckFlow;
        }
        else {
            // Aresta já existe, então é só aumentar o fluxo que passa pela aresta (currentVertex, nextVertex)
            edgeMaps[currentVertex - 1][nextVertex].flow += bottleneckFlow;
        }
    }
}

void WeightedGraph::saveGraphToFile(const std::string& filename, int maxFlow) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << " para escrita." << std::endl;
        return;
    }
    std::cout << "Arquivo aberto com sucesso!" << std::endl;

    outputFile << "Fluxo máximo: " << maxFlow << std::endl;
    // Escrever as arestas e fluxos no arquivo
    for (int currentVertex = 1; currentVertex <= numVertices; ++currentVertex) {
        for (const auto& edge : edgeMaps[currentVertex - 1]) {
            int v = edge.first;
            int capacity = edge.second.capacity;
            int allocatedFlow = edge.second.flow;
            if (allocatedFlow > 0) {
                outputFile << currentVertex << "-->" << v << " / capacidade = " << capacity << ", fluxo = " << allocatedFlow << std::endl;
            }
        }
    }
    outputFile.close();
}

# endif