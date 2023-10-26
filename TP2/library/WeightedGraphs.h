#ifndef WEIGHTEDGRAPHS_H
#define WEIGHTEDGRAPHS_H

#include <bits/stdc++.h> 
using namespace std;
#include <iostream>

class WeightedGraph {
public:

    // Construtor que define o número de vértices
    WeightedGraph(int numVertices);

    // Método para construir o grafo a partir do arquivo lido
    void readGraphFromFile(const std::string& filename);

    // Métodos para adicionar e remover arestas
    void addEdge(int vertexA, int vertexB, float distanceAB);
    void removeEdge(int vertexA, int vertexB);

    // Demais métodos
    int getNumVertices();
    int getNumEdges();
    bool acceptsDijkstra();
    void DijkstraWithVector(int startVertex);
    void DijkstraWithHeap(int startVertex);
    float getMinDistance(int finalVertex);
    std::list<int> minimalPathToVertex(int finalVertex);

private:
    bool hasNegativeDistance = false;
    int numVertices;
    int numEdges;
    int numComponents;
    std::vector<std::list<std::pair<int, float>>> adjacencyList; // Lista de adjacência como vetor dinâmico de listas encadeadas que contém o vértice e a distância
    std::vector<float> dist; // Vetor de distâncias
    std::vector<int> predecessors; // vetor de predecessores
    std::list<int> minimalPath; // lista encadeada para armazenar o menor caminho
};

WeightedGraph::WeightedGraph(int numVertices) {
    this->numVertices = numVertices;
    this->numEdges = 0;
    this->dist.resize(numVertices, std::numeric_limits<float>::infinity()); // Inicializa o vetor de distâncias
    this->predecessors.resize(numVertices, -1);

    adjacencyList.resize(numVertices);
}

void WeightedGraph::addEdge(int vertexA, int vertexB, float distanceAB) {
    adjacencyList[vertexA - 1].push_front(std::make_pair(vertexB, distanceAB));
    numEdges++;
}

void WeightedGraph::removeEdge(int vertexA, int vertexB) {
    for (auto it = adjacencyList[vertexA - 1].begin(); it != adjacencyList[vertexA - 1].end(); ++it) {
        if (it->first == vertexB) {
            adjacencyList[vertexA - 1].erase(it);
            numEdges--;
            break;
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
    
    adjacencyList.resize(numVertices);
    dist.assign(numVertices, std::numeric_limits<float>::infinity()); // Inicializa o vetor de distâncias para cada vértice

    int vertexA, vertexB ;
    float distAB;
    int edgeCount = 0;
    while (inputFile >> vertexA >> vertexB >> distAB) {
        if (distAB < 0) {
            hasNegativeDistance = true;
        }
        addEdge(vertexA, vertexB, distAB);
        addEdge(vertexB, vertexA, distAB);
    }

    inputFile.close();
}

int WeightedGraph::getNumEdges(){
    // Cada aresta está sendo contada 2 vezes (preferi assim pois se fosse um grafo direcionado seria necessário mudar apenas aqui)
    return numEdges/2;
}

int WeightedGraph::getNumVertices(){
    return numVertices;
}

bool WeightedGraph::acceptsDijkstra(){
    return !hasNegativeDistance;
}

void WeightedGraph::DijkstraWithVector(int startVertex) {
    dist.assign(numVertices, std::numeric_limits<float>::infinity());
    predecessors.assign(numVertices, -1); // Inicializa predecessors com -1, indicando que nenhum predecessor foi encontrado ainda
    dist[startVertex - 1] = 0;
    predecessors[startVertex - 1] = 0; // Colocando o predecessor do vértice inicial como 0

    std::set<int> unexploredVertices;
    std::set<int> knownDistances;
    for (int i = 1; i <= numVertices; ++i) {
        unexploredVertices.insert(i);
    }

    // Faz o do primeiro vértice
    for (auto &neighbor : adjacencyList[startVertex - 1]) { //O(grau do vértice)
        int v = neighbor.first;
        float weight = neighbor.second;
        if (dist[v - 1] > dist[startVertex - 1] + weight) {
            dist[v - 1] = dist[startVertex - 1] + weight;
            predecessors[v - 1] = startVertex; // Atualiza o predecessor de v
            knownDistances.insert(v);
            
        }
    }

    // while (!unexploredVertices.empty()) {
    //     int u = -1;
    //     float minDist = std::numeric_limits<float>::infinity();
    //     for (int v : knownDistances) {
    //         if (dist[v - 1] < minDist) {
    //             u = v;
    //             minDist = dist[v - 1];
    //         }
    //     }

    //     if (u == -1) {
    //         break;
    //     }

    //     unexploredVertices.erase(u);
    //     knownDistances.erase(u);

    //     for (auto &neighbor : adjacencyList[u - 1]) {
    //         int v = neighbor.first;
    //         float weight = neighbor.second;

    //         if (dist[v - 1] > dist[u - 1] + weight) {
    //             dist[v - 1] = dist[u - 1] + weight;
    //             predecessors[v - 1] = u;
    //             knownDistances.insert(v);
    //         }
    //     }
    // }



    while (!unexploredVertices.empty()) { // Vai executar n vezes
        int u = *unexploredVertices.begin(); // O(1)
        unexploredVertices.erase(unexploredVertices.begin()); //O(1)
        
        for (auto &neighbor : adjacencyList[u - 1]) { //O(grau do vértice)
            int v = neighbor.first;
            float weight = neighbor.second;

            if (dist[v - 1] > dist[u - 1] + weight) {
                dist[v - 1] = dist[u - 1] + weight;
                predecessors[v - 1] = u; // Atualiza o predecessor de v
                unexploredVertices.insert(v);
            }
        }
    }
}

void WeightedGraph::DijkstraWithHeap(int startVertex) {
    dist.assign(numVertices, std::numeric_limits<float>::infinity());
    predecessors.assign(numVertices, -1); // Inicializa predecessors com -1, indicando que nenhum predecessor foi encontrado ainda
    dist[startVertex - 1] = 0;
    predecessors[startVertex - 1] = 0; // Colocando o predecessor do vértice inicial como 0
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> pq;
    pq.push({0, startVertex});

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        for (auto& neighbor : adjacencyList[u - 1]) {
            int v = neighbor.first;
            float weight = neighbor.second;
            // std::cout << "" << u << 

            if (dist[v - 1] > dist[u - 1] + weight) {
                dist[v - 1] = dist[u - 1] + weight;
                predecessors[v - 1] = u; // Atualiza o predecessor de v
                pq.push({dist[v - 1], v});
            }
        }
    }
}

float WeightedGraph::getMinDistance(int finalVertex) {
    return dist[finalVertex - 1];    
}

std::list<int> WeightedGraph::minimalPathToVertex(int finalVertex){
    minimalPath.clear(); // Limpa a lista encadeada antes de adicionarmos novos elementos
    int currentVertex = finalVertex;

    // Itera até encontrar uma posição vazia em predecessors
    while (currentVertex > 0 && currentVertex <= predecessors.size() && predecessors[currentVertex - 1] > 0) {
        minimalPath.push_front(currentVertex); // Adiciona o vértice atual no início da lista encadeada
        currentVertex = predecessors[currentVertex - 1]; // Atualiza o vértice atual para seu predecessor
    }

    if (currentVertex > 0 && currentVertex <= predecessors.size() && predecessors[currentVertex - 1] == 0) {
        minimalPath.push_front(currentVertex); // Adiciona o último vértice ao caminho mínimo
    }
    return minimalPath;
}

# endif