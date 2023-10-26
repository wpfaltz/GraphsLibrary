#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int graphs[6] = {1,2,3,4,5,6};

int main() {
    // Nome do arquivo de saída
    std::string pathArquivoSaida = "./outputs/Search_Results/Adjacency_List/MeanRunTimes.txt";
    std::ofstream arquivoSaida(pathArquivoSaida);
    for (int i = 0; i < sizeof(graphs) / sizeof(graphs[0]); i++) {
        std::string graphNumberToString = std::to_string(graphs[i]);
        arquivoSaida << ".Grafo " << graphNumberToString << ":" << std::endl;
        
        // Número de linhas a serem lidas
        int numeroDeLinhas = 100;

        // Vetor para armazenar os valores lidos
        std::vector<double> valoresBFS;
        std::string pathToBFSRunTimes = "./outputs/Search_Results/Adjacency_List/Graph_" + graphNumberToString + "/BFS_Elapsed_Times.txt";
        // Abre o arquivo de entrada
        std::ifstream arquivoBFS(pathToBFSRunTimes);
        if (!arquivoBFS.is_open()) {
            std::cerr << "Erro ao abrir o arquivo de entrada." << std::endl;
            return 1;
        }
        // Lê os valores das primeiras 100 linhas
        double valorBFS;
        int linhasLidas = 0;
        while (linhasLidas < numeroDeLinhas && arquivoBFS >> valorBFS) {
            valoresBFS.push_back(valorBFS);
            linhasLidas++;
        }

        // Fecha o arquivo de entrada
        arquivoBFS.close();

        if (valoresBFS.empty()) {
            std::cerr << "Nenhum valor lido do arquivo." << std::endl;
            return 1;
        }

        // Calcula a média dos valores
        double somaBFS = 0;
        for (double v : valoresBFS) {
            somaBFS += v;
        }
        double mediaBFS = somaBFS / valoresBFS.size();

        // Abre o arquivo de saída
        if (!arquivoSaida.is_open()) {
            std::cerr << "Erro ao abrir o arquivo de saída." << std::endl;
            return 1;
        }

        // Escreve a média no arquivo de saída
        arquivoSaida << "   - BFS: " << mediaBFS << " milisegundos"<< std::endl;

        // Fecha o arquivo de saída

        std::cout << "Média de execução da BFS calculada e escrita no arquivo de saída." << std::endl;

        std::vector<double> valoresDFS;
        std::string pathToDFSRunTimes = "./outputs/Search_Results/Adjacency_List/Graph_" + graphNumberToString + "/DFS_Elapsed_Times.txt";
        // Abre o arquivo de entrada
        std::ifstream arquivoDFS(pathToDFSRunTimes);

        if (!arquivoDFS.is_open()) {
            std::cerr << "Erro ao abrir o arquivo de entrada." << std::endl;
            return 1;
        }

        // Lê os valores das primeiras 100 linhas
        double valorDFS;
        linhasLidas = 0;
        while (linhasLidas < numeroDeLinhas && arquivoDFS >> valorDFS) {
            valoresDFS.push_back(valorDFS);
            linhasLidas++;
        }

        // Fecha o arquivo de entrada
        arquivoDFS.close();

        if (valoresDFS.empty()) {
            std::cerr << "Nenhum valor lido do arquivo." << std::endl;
            return 1;
        }

        // Calcula a média dos valores
        double somaDFS = 0;
        for (double v : valoresDFS) {
            somaDFS += v;
        }
        double mediaDFS = somaDFS / valoresDFS.size();

        // Abre o arquivo de saída
        if (!arquivoSaida.is_open()) {
            std::cerr << "Erro ao abrir o arquivo de saída." << std::endl;
            return 1;
        }

        // Escreve a média no arquivo de saída
        arquivoSaida << "   - DFS: " << mediaDFS << " milisegundos"<< std::endl;

        // Fecha o arquivo de saída

        std::cout << "Média de execução da DFS calculada e escrita no arquivo de saída." << std::endl;
    }
    arquivoSaida.close();
    return 0;
}
