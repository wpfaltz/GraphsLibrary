import psutil
import os

memory_stats = psutil.virtual_memory()

def convert_to_GB(amount_in_bytes: int):
    amount_in_GB = amount_in_bytes / (1024 ** 3)
    return amount_in_GB

def convert_to_MB(amount_in_bytes: int):
    print(amount_in_bytes)
    amount_in_MB = amount_in_bytes / (1024 ** 2)
    print(amount_in_MB)
    return amount_in_MB

def tamanho_estimado_matriz_adjacencia(vertices: int, arestas: int):
    tamanho_estimado_bytes = (vertices ** 2) / 8
    tamanho_estimado_gb = convert_to_GB(tamanho_estimado_bytes)
    return tamanho_estimado_gb
    
def tamanho_estimado_lista_adjacencia(vertices: int, arestas: int):
    tamanho_estimado_bytes = 8 * (vertices + (4 * arestas))
    tamanho_estimado_gb = convert_to_GB(tamanho_estimado_bytes)
    return tamanho_estimado_gb

total_memory = memory_stats.total
memory_used = memory_stats.used
free_memory = memory_stats.free

to_convert = [total_memory, memory_used, free_memory]
converted = []

for i in range(len(to_convert)):
    converted.append(convert_to_GB(to_convert[i]))

print("Dados do seu computador:")
print(f". Memória total: {converted[0]:.5f} GB")
print(f"\t- Memória Usada: {converted[1]:.5f} GB")
print(f"\t- Memória Livre: {converted[2]:.5f} GB")

pasta_grafos = r"G:\Outros computadores\Meu laptop\UFRJ\2023.2\Teoria dos Grafos\TP1\graphs"
all_files = os.listdir(pasta_grafos)
print(all_files)
all_files.sort()
only_graphs = []

for i in range(len(all_files)):
    if 'grafo' and 'txt' in all_files[i]:
        only_graphs.append(all_files[i])

print("____________________________________________________________________________")
print("Obs: Todas as estimativas para o tamanho ocupado pelos grafos em cada uma das representações é baseada nas estimativas feitas na lista 2. No caso da matriz de adjacência, foi levado em consideração que cada posição na matriz de adjacência ocupe 1 bit, enquanto que na lista de adjacência cada vértice e cada ponteiro ocupe 8 bytes.")
print("Informações sobre os grafos:")

for i in range(len(only_graphs)):
    path_to_graph = os.path.join(pasta_grafos, only_graphs[i])
    with open(path_to_graph, "r") as fp:
        print(f'. Grafo {i+1}:')
        n = int(fp.readline())
        m = len(fp.readlines())
        print(f'\t - n = {n}')
        print(f'\t - m = {m}')


        tamanho_matriz_adjacencia = tamanho_estimado_matriz_adjacencia(n, m)
        tamanho_lista_adjacencia = tamanho_estimado_lista_adjacencia(n, m)
        print(f'\t - Tamanho estimado a ser ocupado em memória representado como matriz de adjacência: {tamanho_matriz_adjacencia:.3f}GB')
        if tamanho_matriz_adjacencia > converted[2]:
            print('\t\t --> Não será possível representar este grafo por matriz de adjacência neste computador')
        else:
            print('\t\t --> Será possível representar este grafo por matriz de adjacência neste computador sem problemas')
        print(f'\t - Tamanho estimado a ser ocupado em memória representado como lista de adjacência: {tamanho_lista_adjacencia:.3f}GB')
        if tamanho_lista_adjacencia > converted[2]:
            print('\t\t --> Não será possível representar este grafo por lista de adjacência neste computador')
        else:
            print('\t\t --> Será possível representar este grafo por lista de adjacência neste computador sem problemas')