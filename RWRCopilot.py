import networkx as nx
import random as random

class RWRRunner:
    def __init__(self, G, Source) -> None:
        self.G = G
        self.Source = Source

    def RandomWalkWithRestarts(self) -> dict:
        
        self.G = nx.convert_node_labels_to_integers(self.G, first_label=0, ordering='default', label_attribute=None)
        N = self.G.number_of_nodes()
        Nodes = self.G.nodes()
        NodeDict = dict(zip(Nodes, range(N)))
        Source = NodeDict[self.Source]
        RestartProb = 0.15
        P = nx.to_numpy_matrix(self.G)
        P = P / P.sum(axis=0)
        P = (1 - RestartProb) * P + RestartProb * (1 / N)
        A = P.copy()
        for i in range(2, 100):
            A = A * P
            P = P + A
        P = P / P.sum(axis=0)
        return P[:, Source].T.tolist()[0]

    def RandomWalkTest(self) -> dict:
        
        