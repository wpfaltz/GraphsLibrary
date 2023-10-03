#ifndef ADJACENCYMATRIXGRAPH_H
#define ADJACENCYMATRIXGRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <numeric>

class AdjacencyMatrixGraph {
public:
    // Construtor que define o número de vértices e cria a matriz de adjacência
    AdjacencyMatrixGraph(int numVertices);

    // Métodos para adicionar e remover arestas
    void addEdge(int u, int v);
    void removeEdge(int u, int v);

    // Métodos para verificar a existência de arestas e obter informações sobre o grafo
    bool hasEdge(int u, int v);
    int getNumVertices();
    int getNumEdges();
    int getMinDegree();
    int getMaxDegree();
    double getAvgDegree();
    int getMedianDegree();

    // Leitura e escrita do grafo a partir de/para um arquivo
    void readGraphFromFile(const std::string& filename);
    void writeGraphInfoToFile(const std::string& filename);

    // Métodos para busca em largura, busca em profundidade, distâncias e diâmetro
    void BFS(int startVertex, std::vector<int>& distances);
    void BFSNoOutput(int startVertex);
    void BFSWithOutput(int startVertex, const std::string& filename);
    void DFSNoOutput(int startVertex);
    void DFSWithOutput(int startVertex, const std::string& filename);
    int getDistance(int u, int v);
    int getDiameter();

    // Método para encontrar componentes conexas
    std::pair<int, std::vector<std::vector<int>>> findConnectedComponents();

private:
    int numVertices;
    int numEdges;
    std::vector<std::vector<bool>> adjacencyMatrix; // Matriz de adjacência
    std::vector<unsigned char> visited; // Vetor para rastrear os vértices visitados
    std::vector<unsigned char> explored; // Vetor para rastrear vértices explorados (DFS)

    // Função auxiliar para a DFS
    void DFSUtil(int vertex, std::vector<unsigned char>& visited, std::ofstream& outputFile, int depth);

    // Função auxiliar para encontrar componentes conexas
    void DFSConnectedComponents(int vertex, int componentID, std::vector<int>& componentSizes);
    void BFSConnectedComponents(int startVertex, int componentID, std::vector<int>& componentSizes, std::vector<int>& componentVertices);

    // Função auxiliar para calcular a mediana de grau
    int calculateMedianDegree();
};

AdjacencyMatrixGraph::AdjacencyMatrixGraph(int numVertices) {
    this->numVertices = numVertices;
    this->numEdges = 0;

    // Inicialize a matriz de adjacência com falso (sem arestas)
    adjacencyMatrix.assign(numVertices, std::vector<bool>(numVertices, false));
}

void AdjacencyMatrixGraph::addEdge(int u, int v) {
    if (u >= 1 && u <= numVertices && v >= 1 && v <= numVertices && u != v) {
        adjacencyMatrix[u - 1][v - 1] = true;
        adjacencyMatrix[v - 1][u - 1] = true; // Para grafos não direcionados
        numEdges++;
    }
}

void AdjacencyMatrixGraph::removeEdge(int u, int v) {
    if (u >= 1 && u <= numVertices && v >= 1 && v <= numVertices && u != v) {
        adjacencyMatrix[u - 1][v - 1] = false;
        adjacencyMatrix[v - 1][u - 1] = false;
        numEdges--;
    }
}

bool AdjacencyMatrixGraph::hasEdge(int u, int v) {
    if (u >= 1 && u <= numVertices && v >= 1 && v <= numVertices) {
        return adjacencyMatrix[u - 1][v - 1];
    }
    return false;
}

int AdjacencyMatrixGraph::getNumVertices() {
    return numVertices;
}

int AdjacencyMatrixGraph::getNumEdges() {
    return numEdges;
}

int AdjacencyMatrixGraph::getMinDegree() {
    if (numVertices == 0) {
        return 0;
    }

    int minDegree = numVertices;
    for (int u = 0; u < numVertices; u++) {
        int degree = 0;
        for (int v = 0; v < numVertices; v++) {
            if (adjacencyMatrix[u][v]) {
                degree++;
            }
        }
        if (degree < minDegree) {
            minDegree = degree;
        }
    }
    return minDegree;
}

int AdjacencyMatrixGraph::getMaxDegree() {
    if (numVertices == 0) {
        return 0;
    }

    int maxDegree = 0;
    for (int u = 0; u < numVertices; u++) {
        int degree = 0;
        for (int v = 0; v < numVertices; v++) {
            if (adjacencyMatrix[u][v]) {
                degree++;
            }
        }
        if (degree > maxDegree) {
            maxDegree = degree;
        }
    }
    return maxDegree;
}

double AdjacencyMatrixGraph::getAvgDegree() {
    if (numVertices == 0) {
        return 0.0;
    }

    int totalDegree = 0;
    for (int u = 0; u < numVertices; u++) {
        int degree = 0;
        for (int v = 0; v < numVertices; v++) {
            if (adjacencyMatrix[u][v]) {
                degree++;
            }
        }
        totalDegree += degree;
    }

    return static_cast<double>(totalDegree) / numVertices;
}

int AdjacencyMatrixGraph::getMedianDegree() {
    if (numVertices == 0) {
        return 0;
    }

    std::vector<int> degrees(numVertices);
    for (int u = 0; u < numVertices; u++) {
        int degree = 0;
        for (int v = 0; v < numVertices; v++) {
            if (adjacencyMatrix[u][v]) {
                degree++;
            }
        }
        degrees[u] = degree;
    }

    std::sort(degrees.begin(), degrees.end());

    int medianDegree;
    if (numVertices % 2 == 0) {
        medianDegree = (degrees[numVertices / 2 - 1] + degrees[numVertices / 2]) / 2;
    } else {
        medianDegree = degrees[numVertices / 2];
    }

    return medianDegree;
}

void AdjacencyMatrixGraph::readGraphFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    inputFile >> numVertices;
    adjacencyMatrix.assign(numVertices, std::vector<bool>(numVertices, false));

    int u, v;
    int edgeCount = 0;
    while (inputFile >> u >> v) {
        addEdge(u, v);
        edgeCount++;
    }

    numEdges = edgeCount;
    inputFile.close();
}

void AdjacencyMatrixGraph::writeGraphInfoToFile(const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << " para escrita" << std::endl;
        return;
    }

    outputFile << "Número de vértices: " << numVertices << std::endl;
    outputFile << "Número de arestas: " << numEdges << std::endl;

    if (numVertices > 0) {
        int minDegree = getMinDegree();
        int maxDegree = getMaxDegree();
        double avgDegree = getAvgDegree();
        int medianDegree = getMedianDegree();
        
        outputFile << "Grau Mínimo: " << minDegree << std::endl;
        outputFile << "Grau Máximo: " << maxDegree << std::endl;
        outputFile << "Grau Médio: " << avgDegree << std::endl;
        outputFile << "Mediana de Grau: " << medianDegree << std::endl;

        std::pair<int, std::vector<std::vector<int>>> componentInfo = findConnectedComponents();
        int numComponents = componentInfo.first;
        std::vector<std::vector<int>> componentList = componentInfo.second;

        outputFile << "Número de Componentes Conexas: " << numComponents << std::endl;
        outputFile << "Tamanho das Componentes Conexas:" << std::endl;
        for (int i = 0; i < numComponents; ++i) {
            outputFile << " - Componente " << i + 1 << ": " << componentList[i].size() << " vértices" << std::endl;
        }

        outputFile << "Lista de Vértices por Componente Conexa:" << std::endl;
        for (int i = 0; i < numComponents; ++i) {
            outputFile << " - Componente " << i + 1 << ": ";
            for (int vertex : componentList[i]) {
                outputFile << vertex << " ";
            }
            outputFile << std::endl;
        }
    }

    outputFile.close();
}

void AdjacencyMatrixGraph::BFS(int startVertex, std::vector<int>& distances) {
    visited.assign(numVertices, 0);
    std::queue<int> queue;

    visited[startVertex - 1] = 1;
    queue.push(startVertex);
    distances[startVertex - 1] = 0;

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
            if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                queue.push(neighbor);
                distances[neighbor - 1] = distances[currentVertex - 1] + 1;
            }
        }
    }
}

void AdjacencyMatrixGraph::BFSNoOutput(int startVertex) {
    visited.assign(numVertices, 0);
    std::queue<int> queue;
    std::vector<int> levels(numVertices, -1);

    levels[startVertex - 1] = 0;
    visited[startVertex - 1] = 1;
    queue.push(startVertex);

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
            if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                queue.push(neighbor);
            }
        }
    }
}

void AdjacencyMatrixGraph::BFSWithOutput(int startVertex, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    visited.assign(numVertices, 0);
    std::queue<int> queue;
    std::vector<int> levels(numVertices, -1);

    levels[startVertex - 1] = 0;
    visited[startVertex - 1] = 1;
    queue.push(startVertex);
    outputFile << "Vértice: " << startVertex << ", Nível: " << levels[startVertex - 1] << ", Nó raiz" << std::endl;

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
            if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                queue.push(neighbor);
                outputFile << "Vértice: " << neighbor << ", Nível: " << levels[neighbor - 1] << ", Pai: " << currentVertex << std::endl;
            }
        }
    }

    outputFile.close();
}

void AdjacencyMatrixGraph::DFSNoOutput(int startVertex) {
    explored.assign(numVertices, 0);
    std::stack<int> exploreStack;
    std::vector<int> fathers(numVertices, 0);
    std::vector<int> levels(numVertices, -1);

    levels[startVertex - 1] = 0;
    fathers[startVertex - 1] = 0;
    exploreStack.push(startVertex);

    while (!exploreStack.empty()) {
        int currentVertex = exploreStack.top();
        exploreStack.pop();

        if (explored[currentVertex - 1] == 0) {
            explored[currentVertex - 1] = 1;
            for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
                if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !explored[neighbor - 1]) {
                    levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                    fathers[neighbor - 1] = currentVertex;
                    exploreStack.push(neighbor);
                }
            }
        }
    }
}

void AdjacencyMatrixGraph::DFSWithOutput(int startVertex, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    explored.assign(numVertices, 0);
    std::stack<int> exploreStack;
    std::vector<int> fathers(numVertices, 0);
    std::vector<int> levels(numVertices, -1);

    levels[startVertex - 1] = 0;
    fathers[startVertex - 1] = 0;
    exploreStack.push(startVertex);

    while (!exploreStack.empty()) {
        int currentVertex = exploreStack.top();
        exploreStack.pop();

        if (fathers[currentVertex - 1] == 0 && explored[currentVertex - 1] == 0) {
            outputFile << "Vértice: " << currentVertex << ", Nível: " << levels[currentVertex - 1] << ", Nó raiz" << std::endl;
        }
        else if (explored[currentVertex - 1] == 0) {
            outputFile << "Vértice: " << currentVertex << ", Nível: " << levels[currentVertex - 1] << ", Pai: " << fathers[currentVertex - 1] << std::endl;
        }

        if (explored[currentVertex - 1] == 0) {
            explored[currentVertex - 1] = 1;
            for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
                if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !explored[neighbor - 1]) {
                    levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                    fathers[neighbor - 1] = currentVertex;
                    exploreStack.push(neighbor);
                }
            }
        }
    }

    outputFile.close();
}

int AdjacencyMatrixGraph::getDistance(int startVertex, int finalVertex) {
    if (startVertex < 1 || startVertex > numVertices || finalVertex < 1 || finalVertex > numVertices) {
        std::cerr << "Vértices inválidos." << std::endl;
        return -1;
    }

    visited.assign(numVertices, 0);
    std::vector<int> distance(numVertices, -1);
    std::queue<int> bfsQueue;

    visited[startVertex - 1] = 1;
    bfsQueue.push(startVertex);
    distance[startVertex - 1] = 0;

    while (!bfsQueue.empty()) {
        int currentVertex = bfsQueue.front();
        bfsQueue.pop();

        if (currentVertex == finalVertex) {
            return distance[finalVertex - 1];
        }

        for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
            if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                bfsQueue.push(neighbor);
                distance[neighbor - 1] = distance[currentVertex - 1] + 1;
            }
        }
    }

    return -1;
}

int AdjacencyMatrixGraph::getDiameter() {
    int diameter = 0;

    std::vector<std::thread> threads;
    std::vector<int> threadResults(numVertices, -1);
    std::mutex mutex;

    for (int startVertex = 1; startVertex <= numVertices; ++startVertex) {
        threads.emplace_back([startVertex, &threadResults, &mutex, this]() {
            std::vector<int> distances(numVertices, -1);
            BFS(startVertex, distances);

            int highestShortestPath = *std::max_element(distances.begin(), distances.end());

            std::lock_guard<std::mutex> lock(mutex);
            if (highestShortestPath > threadResults[startVertex - 1]) {
                threadResults[startVertex - 1] = highestShortestPath;
            }
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    diameter = *std::max_element(threadResults.begin(), threadResults.end());

    return diameter;
}

std::pair<int, std::vector<std::vector<int>>> AdjacencyMatrixGraph::findConnectedComponents() {
    int componentID = 0;
    std::vector<int> componentSizes;
    std::vector<std::vector<int>> componentList;
    visited.assign(numVertices, 0);

    for (int vertex = 1; vertex <= numVertices; vertex++) {
        if (!visited[vertex - 1]) {
            componentID++;
            componentSizes.push_back(0);
            componentList.push_back(std::vector<int>());

            if (componentID > 1) {
                BFSConnectedComponents(vertex, componentID, componentSizes, componentList.back());
            } else {
                DFSConnectedComponents(vertex, componentID, componentSizes);
            }
        }
    }

    return std::make_pair(componentID, componentList);
}

void AdjacencyMatrixGraph::DFSConnectedComponents(int vertex, int componentID, std::vector<int>& componentSizes) {
    visited[vertex - 1] = 1;
    componentSizes[componentID - 1]++;
    for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
        if (adjacencyMatrix[vertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
            DFSConnectedComponents(neighbor, componentID, componentSizes);
        }
    }
}

void AdjacencyMatrixGraph::BFSConnectedComponents(int startVertex, int componentID, std::vector<int>& componentSizes, std::vector<int>& componentVertices) {
    visited[startVertex - 1] = 1;
    std::queue<int> bfsQueue;
    bfsQueue.push(startVertex);
    componentSizes[componentID - 1]++;
    componentVertices.push_back(startVertex);

    while (!bfsQueue.empty()) {
        int currentVertex = bfsQueue.front();
        bfsQueue.pop();

        for (int neighbor = 1; neighbor <= numVertices; neighbor++) {
            if (adjacencyMatrix[currentVertex - 1][neighbor - 1] && !visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                bfsQueue.push(neighbor);
                componentSizes[componentID - 1]++;
                componentVertices.push_back(neighbor);
            }
        }
    }
}

int AdjacencyMatrixGraph::calculateMedianDegree() {
    if (numVertices == 0) {
        return 0;
    }

    std::vector<int> degrees(numVertices);
    for (int u = 0; u < numVertices; u++) {
        int degree = 0;
        for (int v = 0; v < numVertices; v++) {
            if (adjacencyMatrix[u][v]) {
                degree++;
            }
        }
        degrees[u] = degree;
    }

    std::sort(degrees.begin(), degrees.end());

    int medianDegree;
    if (numVertices % 2 == 0) {
        medianDegree = (degrees[numVertices / 2 - 1] + degrees[numVertices / 2]) / 2;
    } else {
        medianDegree = degrees[numVertices / 2];
    }

    return medianDegree;
}

#endif  // ADJACENCYMATRIXGRAPH_H

#ifndef ADJACENCYLISTGRAPH_H
#define ADJACENCYLISTGRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <thread>
#include <mutex>
#include <string>
#include <algorithm>
#include <numeric>

class AdjacencyListGraph {
public:
    // Construtor que define o número de vértices e cria a lista de adjacência
    AdjacencyListGraph(int numVertices);

    // Métodos para adicionar e remover arestas
    void addEdge(int u, int v);
    void removeEdge(int u, int v);

    // Métodos para verificar a existência de arestas e obter informações sobre o grafo
    bool hasEdge(int u, int v);
    int getNumVertices();
    int getNumEdges();
    int getMinDegree();
    int getMaxDegree();
    double getAvgDegree();
    int getMedianDegree();

    // Leitura e escrita do grafo a partir de/para um arquivo
    void readGraphFromFile(const std::string& filename);
    void writeGraphInfoToFile(const std::string& filename);

    // Métodos para busca em largura, busca em profundidade, distâncias e diâmetro
    int BFS(int startVertex);
    void BFSNoOutput(int startVertex);
    void BFSWithOutput(int startVertex, const std::string& filename);
    void DFSNoOutput(int startVertex);
    void DFSWithOutput(int startVertex, const std::string& filename);
    int getDistance(int u, int v);
    int getDiameter();

    // Método para encontrar componentes conexas
    std::pair<int, std::vector<std::vector<int>>> findConnectedComponents();

private:
    int numVertices;
    int numEdges;
    int numComponents;
    std::vector<std::list<int>> adjacencyList; // Lista de adjacência como vetor dinâmico de listas encadeadas
    std::vector<unsigned char> visited; // Vetor para rastrear os vértices visitados
    std::vector<unsigned char> explored; // Vetor para rastrear vértices explorados (DFS)

    // Função auxiliar para a DFS
    void DFSUtil(int vertex, std::vector<unsigned char>& visited, std::ofstream& outputFile, int depth);

    // Função auxiliar para encontrar componentes conexas
    void DFSConnectedComponents(int vertex, int componentID, std::vector<int>& componentSizes);
    void BFSConnectedComponents(int startVertex, int componentID, std::vector<int>& componentSizes, std::vector<int>& componentVertices);

    // Função auxiliar para calcular a mediana de grau
    int calculateMedianDegree();
};

AdjacencyListGraph::AdjacencyListGraph(int numVertices) {
    this->numVertices = numVertices;
    this->numEdges = 0;

    // Inicialize a lista de adjacência com um vetor de listas vazias
    adjacencyList.resize(numVertices);
}

void AdjacencyListGraph::addEdge(int u, int v) { // OK
    if (u >= 0 && u < numVertices && v >= 0 && v < numVertices && u != v) {
        adjacencyList[u].push_back(v);
        adjacencyList[v].push_back(u); // Para grafos não direcionados
        numEdges++;
    }
}

void AdjacencyListGraph::removeEdge(int u, int v) {
    if (u >= 0 && u < numVertices && v >= 0 && v < numVertices && u != v) {
        // Remova v da lista de adjacência de u e vice-versa
        adjacencyList[u - 1].remove(v);
        adjacencyList[v - 1].remove(u);
        numEdges--;
    }
}

bool AdjacencyListGraph::hasEdge(int u, int v) { // OK
    if (u >= 0 && u < numVertices && v>= 0 && v < numVertices) {
        // Vertifique se 'v' está na lista de adjacência de 'u'
        for (int neighbor : adjacencyList[u - 1]) {
            if (neighbor == v) {
                return true; // Aresta encontrada
            }
        }
    }
    return false; // Aresta não encontrada ou vértices fora dos limites
}

int AdjacencyListGraph::getNumVertices() { // OK
    return numVertices;
}

int AdjacencyListGraph::getNumEdges() { // OK
    return numEdges;
}

int AdjacencyListGraph::getMinDegree(){ //OK
    if (numVertices == 0) {
        return 0; // Retorna 0 se não houver vértices no grafo
    }

    int minDegree = numVertices; // Inicialize com um valor alto

    // Itere sobre todos os vértices e atualize o valor mínimo
    for (int u = 0; u < numVertices; u++) {
        int degree = adjacencyList[u].size();
        if (degree < minDegree) {
            minDegree = degree;
        }
    }

    return minDegree;
}

int AdjacencyListGraph::getMaxDegree() { // OK
    if (numVertices == 0) {
        return 0; // Retorna 0 se não houver vértices no grafo
    }

    int maxDegree = 0;

    // Itere sobre todos os vértices e atualize o valor máximo
    for (int u = 0; u < numVertices; u++) {
        int degree = adjacencyList[u].size();
        if (degree > maxDegree) {
            maxDegree = degree;
        }
    }

    return maxDegree;
}

double AdjacencyListGraph::getAvgDegree() { // OK
    if (numVertices == 0) {
        return 0.0; // Retorna 0.0 se não houver vértices no grafo
    }

    int totalDegree = 0;

    // Itere sobre todos os vértices e some os graus
    for (int u = 0; u < numVertices; u++) {
        totalDegree += adjacencyList[u].size();
    }

    return static_cast<double>(totalDegree) / numVertices;
}

int AdjacencyListGraph::getMedianDegree() { // OK
    if (numVertices == 0) {
        return 0; // Retorna 0 se não houver vértices no grafo
    }

    // Crie um vetor para armazenar os graus
    std::vector<int> degrees(numVertices);

    // Preencha o vetor com os graus dos vértices
    for (int u = 0; u < numVertices; u++) {
        degrees[u] = adjacencyList[u].size();
    }

    // Ordene o vetor de graus
    std::sort(degrees.begin(), degrees.end());

    // Calcule a mediana
    int medianDegree;
    if (numVertices % 2 == 0) {
        medianDegree = (degrees[numVertices / 2 - 1] + degrees[numVertices / 2]) / 2;
    } else {
        medianDegree = degrees[numVertices / 2];
    }

    return medianDegree;
}

void AdjacencyListGraph::readGraphFromFile(const std::string& filename) { // OK
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    inputFile >> numVertices;
    // std::cerr << "Número de vértices: " << numVertices << std::endl;

    // Redimensione a lista de adjacência com base no número de vértices
    adjacencyList.resize(numVertices);

    int u, v;
    int edgeCount = 0;
    while (inputFile >> u >> v) {
        // Adicione a aresta à lista de adjacência (grafo não direcionado)
        // std::cerr << "Aresta " << edgeCount+1 << " lida: " << u << " --> " << v << std::endl;
        adjacencyList[u-1].insert(adjacencyList[u-1].begin(), v);
        // std::cerr << "Vértice " << v << " adicionado à lista encadeada do vértice " << u << std::endl;
        adjacencyList[v-1].insert(adjacencyList[v-1].begin(), u);
        // std::cerr << "Vértice " << u << " adicionado à lista encadeada do vértice " << v << std::endl;
        edgeCount++;
    }

    numEdges = edgeCount;
    // std::cerr << "Total de arestas: " << numEdges << std::endl;
    inputFile.close();
    // std::cerr << "Arquivo fechado" << std::endl;
}

void AdjacencyListGraph::writeGraphInfoToFile(const std::string& filename) { // OK
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << " para escrita" << std::endl;
        return;
    }

    // Escreva informações sobre o grafo no arquivo
    outputFile << "Número de vértices: " << numVertices << std::endl;
    outputFile << "Número de arestas: " << numEdges << std::endl;

    if (numVertices > 0) {
        // Calcule e escreva grau mínimo, máximo, médio e mediana de grau
        int minDegree = getMinDegree();
        int maxDegree = getMaxDegree();
        double avgDegree = getAvgDegree();
        int medianDegree = getMedianDegree();
        
        outputFile << "Grau Mínimo: " << minDegree << std::endl;
        outputFile << "Grau Máximo: " << maxDegree << std::endl;
        outputFile << "Grau Médio: " << avgDegree << std::endl;
        outputFile << "Mediana de Grau: " << medianDegree << std::endl;

        // Encontre informações sobre componentes conexas
        std::pair<int, std::vector<std::vector<int>>> componentInfo = findConnectedComponents();
        int numComponents = componentInfo.first;
        std::vector<std::vector<int>> componentList = componentInfo.second;

        outputFile << "Número de Componentes Conexas: " << numComponents << std::endl;
        outputFile << "Tamanho das Componentes Conexas:" << std::endl;
        // Escreva os tamanhos das componentes
        for (int i = 0; i < numComponents; ++i) {
            outputFile << " - Componente " << i + 1 << ": " << componentList[i].size() << " vértices" << std::endl;
        }

        outputFile << "Lista de Vértices por Componente Conexa:" << std::endl;
        // Escreva os vértices por componente
        for (int i = 0; i < numComponents; ++i) {
            outputFile << " - Componente " << i + 1 << ": ";
            for (int vertex : componentList[i]) {
                outputFile << vertex << " ";
            }
            outputFile << std::endl;
        }

        // int diameter = getDiameter();
        // outputFile << "Diâmetro: " << diameter << std::endl;
    }

    outputFile.close();
}

int AdjacencyListGraph::BFS(int startVertex) { // OK
    visited.assign(numVertices, 0);
    std::queue<int> queue;

    std::vector<int> levels;
    levels.assign(numVertices, -1);

    levels[startVertex - 1] = 0;

    visited[startVertex - 1] = 1;
    queue.push(startVertex);

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        // Processar os vértices adjacentes não visitados
        for (int neighbor : adjacencyList[currentVertex - 1]) {
            if (visited[neighbor - 1] == 0) {
                visited[neighbor - 1] = 1;
                levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                queue.push(neighbor);
            }
        }
    }
    auto maior_elemento = std::max_element(levels.begin(), levels.end());

    if (maior_elemento != levels.end()) {
        int maiorDistancia = *maior_elemento;
        return maiorDistancia;
    } else {
        std::cout << "O vetor está vazio." << std::endl;
        return -1;
    }
}

void AdjacencyListGraph::BFSNoOutput(int startVertex) {
    visited.assign(numVertices, 0);
    std::queue<int> queue;

    std::vector<int> levels;
    levels.assign(numVertices, -1);

    // Inicialize o nível do vértice de partida como zero
    levels[startVertex - 1] = 0;

    visited[startVertex - 1] = 1;
    queue.push(startVertex);

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        // Processar os vértices adjacentes não visitados
        for (int neighbor : adjacencyList[currentVertex - 1]) {
            if (visited[neighbor - 1] == 0) {
                visited[neighbor - 1] = 1;
                levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                queue.push(neighbor);
            }
        }
    }
    auto maior_elemento = std::max_element(levels.begin(), levels.end());

    if (maior_elemento != levels.end()) {
        int maiorDistancia = *maior_elemento;
        std::cout << "O maior valor de distäncia é " << maiorDistancia << std::endl;
    } else {
        std::cout << "O vetor está vazio." << std::endl;
    }
}

void AdjacencyListGraph::BFSWithOutput(int startVertex, const std::string& filename) { // OK
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }

    visited.assign(numVertices, 0);
    std::queue<int> queue;

    std::vector<int> levels;
    levels.assign(numVertices, -1);

    // Inicialize o nível do vértice de partida como zero
    levels[startVertex - 1] = 0;

    visited[startVertex - 1] = 1;
    queue.push(startVertex);
    outputFile << "Vértice: " << startVertex << ", Nível: " << levels[startVertex - 1] << ", Nó raiz" << std::endl;

    while (!queue.empty()) {
        int currentVertex = queue.front();
        queue.pop();

        // Processar os vértices adjacentes não visitados
        for (int neighbor : adjacencyList[currentVertex - 1]) {
            if (visited[neighbor - 1] == 0) {
                visited[neighbor - 1] = 1;
                levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                queue.push(neighbor);
                outputFile << "Vértice: " << neighbor << ", Nível: " << levels[neighbor - 1] << ", Pai: " << currentVertex << std::endl;
            }
        }
    }

    outputFile.close();
}

void AdjacencyListGraph::DFSNoOutput(int startVertex) {
    explored.assign(numVertices, 0);
    std::stack<int> exploreStack;
    std::vector<int> fathers;
    fathers.assign(numVertices, 0);
    std::vector<int> levels;
    levels.assign(numVertices, -1);

    // Inicialize o nível do vértice de partida como 0
    levels[startVertex - 1] = 0;
    exploreStack.push(startVertex);

    while (!exploreStack.empty()) {
        int currentVertex = exploreStack.top();
        exploreStack.pop();        

        if (explored[currentVertex - 1] == 0) {
            explored[currentVertex - 1] = 1;
            for (int neighbor : adjacencyList[currentVertex - 1]) {
                    levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                    fathers[neighbor - 1] = currentVertex;
                    exploreStack.push(neighbor);                 
                }
        }
    }
}

void AdjacencyListGraph::DFSWithOutput(int startVertex, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return;
    }
    explored.assign(numVertices, 0);
    std::stack<int> exploreStack;
    std::vector<int> fathers;
    fathers.assign(numVertices, 0);
    std::vector<int> levels;
    levels.assign(numVertices, -1);

    // Inicialize o nível do vértice de partida como 0
    levels[startVertex - 1] = 0;
    exploreStack.push(startVertex);

    while (!exploreStack.empty()) {
        int currentVertex = exploreStack.top();
        exploreStack.pop();

        if (fathers[currentVertex - 1] == 0 && explored[currentVertex - 1] == 0) {
            outputFile << "Vértice: " << currentVertex << ", Nível: " << levels[currentVertex - 1] << ", Nó raiz" << std::endl;
        }
        else if (explored[currentVertex - 1] == 0){
            outputFile << "Vértice: " << currentVertex << ", Nível: " << levels[currentVertex - 1] << ", Pai: " << fathers[currentVertex - 1] << std::endl;
        }
        

        if (explored[currentVertex - 1] == 0) {
            explored[currentVertex - 1] = 1;
            for (int neighbor : adjacencyList[currentVertex - 1]) {
                    levels[neighbor - 1] = levels[currentVertex - 1] + 1;
                    fathers[neighbor - 1] = currentVertex;
                    exploreStack.push(neighbor);                 
                }
        }
    }
    outputFile.close(); // Feche o arquivo de saída
}

int AdjacencyListGraph::getDistance(int startVertex, int finalVertex) { // OK
    // Verifique se os vértices são válidos
    if (startVertex < 1 || startVertex > numVertices || finalVertex < 1 || finalVertex > numVertices) {
        std::cerr << "Vértices inválidos." << std::endl;
        return -1; // Valor de retorno inválido para indicar erro
    }

    visited.assign(numVertices, 0);

    // Vetor de distâncias a partir do vértice 'u'
    std::vector<int> distance(numVertices, -1);

    // Crie uma fila para a BFS
    std::queue<int> bfsQueue;
    
    visited[startVertex - 1] = 1; 
    bfsQueue.push(startVertex);
    distance[startVertex - 1] = 0;


    while (!bfsQueue.empty()) {
        int currentVertex = bfsQueue.front();
        bfsQueue.pop();

        // Verifique se alcançamos o vértice 'v'
        if (currentVertex == finalVertex) {
            return distance[finalVertex - 1];
        }

        // Explore todos os vértices adjacentes
        for (int neighbor : adjacencyList[currentVertex - 1]) {
            if (!visited[neighbor - 1]) {
                visited[neighbor - 1] = 1;
                // Se o vértice não foi visitado, atualize a distância e coloque-o na fila
                bfsQueue.push(neighbor);
                distance[neighbor - 1] = distance[currentVertex - 1] + 1;
            }
        }
    }

    // Se não foi possível alcançar o vértice 'v' a partir de 'u', retorne -1
    return -1;
}

int AdjacencyListGraph::getDiameter() { // OK

    // Inicialize a variável para armazenar o diâmetro
    int diameter = 0;

    for (int startVertex = 1; startVertex <= numVertices; ++startVertex) {
            int maxLevel = BFS(startVertex);
            if (maxLevel > diameter) {
                diameter = maxLevel;
            }
    }
    return diameter;
}

std::pair<int, std::vector<std::vector<int>>> AdjacencyListGraph::findConnectedComponents() { // OK
    std::vector<std::vector<int>> components;
    std::vector<int> componentSizes(numVertices, 0);
    int componentID = 0;

    visited.assign(numVertices, 0); // Reinicia o vetor de visitados

    for (int vertex = 1; vertex <= numVertices; ++vertex) {
        if (visited[vertex - 1] == 0) {
            // Encontrou um novo componente conexo
            components.emplace_back(); // Adiciona um novo componente
            BFSConnectedComponents(vertex, componentID, componentSizes, components.back());
            componentID++;
        }
    }

    return std::make_pair(componentID, components);
}

void AdjacencyListGraph::DFSUtil(int vertex, std::vector<unsigned char>& visited, std::ofstream& outputFile, int depth) { 
    // Marque o vértice como visitado
    visited[vertex - 1] = true;

    // Escreva informações sobre o vértice no arquivo de saída, incluindo a profundidade
    outputFile << "Vértice " << vertex+1 << " visitado com profundidade " << depth << std::endl;

    // Recorra todos os vértices adjacentes não visitados
    for (const int& neighbor : adjacencyList[vertex]) {
        if (!visited[neighbor - 1]) {
            // Chame a DFS recursivamente para o vértice adjacente com profundidade + 1
            DFSUtil(neighbor, visited, outputFile, depth + 1);
        }
    }
}

void AdjacencyListGraph::BFSConnectedComponents(int startVertex, int componentID, std::vector<int>& componentSizes, std::vector<int>& componentVertices) { //OK
    visited[startVertex - 1] = 1; // Marca o vértice inicial como visitado
    componentSizes[componentID]++;
    componentVertices.push_back(startVertex);

    std::queue<int> bfsQueue;
    bfsQueue.push(startVertex);

    while (!bfsQueue.empty()) {
        int currentVertex = bfsQueue.front();
        bfsQueue.pop();

        for (int neighbor : adjacencyList[currentVertex - 1]) {
            if (visited[neighbor - 1] == 0) {
                visited[neighbor - 1] = 1;
                componentSizes[componentID]++;
                componentVertices.push_back(neighbor);
                bfsQueue.push(neighbor);
            }
        }
    }
}

#endif