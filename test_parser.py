from parser_module import *

#default threshold is 1%
#make sure these files are in the expected directory
data = Annotation_Parser('go.obo', 'goa_human.gaf')
GO_graph = data.GO_graph
print('Number of GO terms:', GO_graph.number_of_nodes())

annotations = data.annotations

ancestors = data.ancestors_by_GO

print('\nThreshold = 1%')
S = data.S
S_length = {key: len(value) for key, value in S.items()}
print('Terms in S:' , S_length)

T = data.T
T_length = {key: len(value) for key, value in T.items()}
print('Terms in T:', T_length)

#data.stacked_plot()

#Changing the threshold
data.threshold = 10

print('\nThreshold = 10%')
S_10 = data.S
S_length_10 = {key: len(value) for key, value in S_10.items()}
print('Terms in S:' , S_length_10)

T_10 = data.T
T_length_10 = {key: len(value) for key, value in T_10.items()}
print('Terms in T:', T_length_10)

#data.stacked_plot()