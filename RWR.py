import networkx as nx
import numpy as np

G = nx.karate_club_graph()

adj = nx.adjacency_matrix(G, None, None, None)

Sources = [1, 12, 17, 33]

s = np.zeros(adj.shape[0])

for Source in Sources:
    s[Source] = 1 / len(Sources)

diag = np.diag(np.sum(adj, axis=1).A1)

inverseDiag = np.linalg.inv(diag)

prod = inverseDiag * adj

prodTranspose = prod.transpose()

identityMatrix = np.identity(adj.shape[0])

vectorStationaryProbabilities = np.linalg.inv(identityMatrix - ((1 - 0.15) * prodTranspose))

a = np.array(vectorStationaryProbabilities)

b = np.array(s)

vectorStationaryProbabilities = np.dot(a, b)

print(vectorStationaryProbabilities)

dictVectorStationaryProbabilities = dict(zip(G.nodes(), vectorStationaryProbabilities))

value_key_pairs = ((value, key) for (key,value) in dictVectorStationaryProbabilities.items())
sorted_value_key_pairs = sorted(value_key_pairs, reverse=True)

print(sorted_value_key_pairs)