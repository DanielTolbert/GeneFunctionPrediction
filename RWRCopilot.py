import networkx as nx
import random as random

class RWRRunner:
    def __init__(self, G: nx.Graph(), Source) -> None:
        self.G = G
        self.Source = Source

    def RandomWalkWithRestarts(self) -> dict:
        
        self.G = nx.convert_node_labels_to_integers(self.G, first_label=0, ordering='default', label_attribute=None)
        N = self.G.number_of_nodes()
        Nodes = self.G.nodes()
        # NodeDict = dict(zip(Nodes, range(N)))
        # for node, n in NodeDict.items():
            # print("{}   ({})".format(node, n))
        # Source = NodeDict[self.Source]
        Source = self.Source
        RestartProb = 0.15
        P = nx.to_numpy_matrix(self.G)
        # print(P.sum(axis=0))
        P = P / P.sum(axis=0)
        P = (1 - RestartProb) * P + RestartProb * (1 / N)
        A = P.copy()
        for i in range(2, 100):
            A = A * P
            P = P + A
        P = P / P.sum(axis=0)
        return P[:, Source].T.tolist()[0]


OldG = nx.karate_club_graph()

# f1 = csv.reader(open("WHPINEdges.csv"))


# with open("EGFR1Edges.csv") as file:
#     edges = file.readlines()
#     for edge in edges:
#         OldG.add_edge(edge.split(",")[0], edge.split(",")[1])
# print(OldG.number_of_nodes())
# print(OldG.number_of_edges())

# for node in OldG.nodes():
    # print(node)

runner = RWRRunner(OldG, {34, 22, 32})
print(runner.RandomWalkWithRestarts())


