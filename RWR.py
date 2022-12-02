import networkx as nx
import numpy as np
from parser_module import *
import time
# To run, python .\RWR.py > your_output_file_name.txt

t = time.time()
data = Annotation_Parser('go.obo', 'goa_human.gaf', 'human_interactome_weighted.txt', threshold=1)
G = data.interactome

# The first run will take about 14 minutes, but once the run is complete, will make subsequent runs faster.
# Uncomment lines 13 to 29 for the first run, and comment out line 33 for the first run.
# adj = nx.adjacency_matrix(G, None, None, None)

# s = np.zeros(adj.shape[0])

# diag = np.diag(np.sum(adj, axis=1).A1)

# inverseDiag = np.linalg.pinv(diag)

# prod = inverseDiag * adj

# prodTranspose = prod.transpose()

# identityMatrix = np.identity(adj.shape[0])

# vectorStationaryProbabilities = np.linalg.inv(identityMatrix - ((1 - 0.5) * prodTranspose))

# np.save('vectorStationaryProbabilities.npy', vectorStationaryProbabilities)

# After the first run, re comment lines 13 to 29, and uncomment line 33

vectorStationaryProbabilities = np.load('vectorStationaryProbabilities.npy')
# For each go term, change the value below, I suggest using a different output file each time.
pos, neg = data.examples('GO:0048523')
# print(pos) # In order to randomly set some positive examples to 0, set newVar equal to the new value of pos with some positive values set to 0.
# newVar = 
# Replace pos in line 38 with newVar
s = pos / np.sum(pos)


a = np.array(vectorStationaryProbabilities)

b = np.array(s)

vectorStationaryProbabilities = np.dot(a, b) * .5

dictVectorStationaryProbabilities = dict(zip(G.nodes(), vectorStationaryProbabilities))

value_key_pairs = ((value, key) for (key,value) in dictVectorStationaryProbabilities.items())
sorted_value_key_pairs = sorted(value_key_pairs, reverse=True)
t2 = time.time()
print("Run took seconds: ", t2 - t)
print("Below is the positive examples: ")
nodeDict = dict(zip(data.genes, range(len(G.nodes()))))
for i in range(len(sorted_value_key_pairs)):
    if pos[nodeDict.get(sorted_value_key_pairs[i][1])] == 1:
        print("Node: ", sorted_value_key_pairs[i][1], " Score: ", sorted_value_key_pairs[i][0], " Positive Example")
    elif neg[nodeDict.get(sorted_value_key_pairs[i][1])] == 1:
        print("Node: ", sorted_value_key_pairs[i][1], " Score: ", sorted_value_key_pairs[i][0], " Negative Example")
    else:
        print("Node: ", sorted_value_key_pairs[i][1], " Score: ", sorted_value_key_pairs[i][0], " Non Example")
print(np.sum(vectorStationaryProbabilities))