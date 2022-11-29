from parser_module import *
import numpy as np

#default threshold is 1%
#make sure these files are in the expected directory
data = Annotation_Parser('go.obo', 'goa_human.gaf', 'human_interactome_weighted.txt', threshold=10)
GO_graph = data.GO_graph
genes = data.genes
print('Number of genes in interactome:', len(genes))

T_10 = data.T
T_length_10 = {key: len(value) for key, value in T_10.items()}
print('Terms in T:', T_length_10)
print(T_10)
#data.stacked_plot()

for namespace in ROOT_TERMS:
    for term in T_10[namespace]:
        pos, neg = data.examples(term)
        print(term, np.sum(pos), np.sum(neg))