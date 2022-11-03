import networkx as nx
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv

neighbors = {}
node_index = {}

def human_ppi_parser():
    ppi_file = open("pathlinker-human-network.txt",'r')
    ppi_lines = ppi_file.readlines()
    g = nx.DiGraph()

    skip = 0
    index = 0
    for lines in ppi_lines:
        if skip == 0:
            skip += 1
        else:
            data = lines.split('\t')
            if data[0] not in g.nodes():
                g.add_node(data[0])
                neighbors[data[0]] = []
                node_index[data[0]] = index
                index += 1
            if data[1] not in g.nodes():
                g.add_node(data[1])
                neighbors[data[1]] = []
                node_index[data[1]] = index
                index += 1
            
            g.add_edge(data[0], data[1])
            neighbors[data[0]].append(data[1])

    return g


#given a graph and alpha: return a reguralized lapacian
def rl(g,alpha):
    #W : adjacency matrix of G
    node_len = len(list(g.nodes))
    row = []
    coln = []
    data = []
    for edge in g.edges():
        row.append(node_index[edge[0]])
        coln.append(node_index[edge[1]])
        data.append(1)
    sparse_w = csc_matrix((data, (row, coln)),shape=(node_len,node_len))
    #D: diagonal matrix where Duu = Sum[v] w_uv , for every node u in G
    row = []
    coln = []
    data = []
    for u in neighbors:
        if len(neighbors[u]) > 0:
            row.append(node_index[u])
            coln.append(node_index[u])
            data.append(len(neighbors[u]))
    sparse_d = csc_matrix((data, (row, coln)),shape=(node_len,node_len))
    #W˜ = D^(−1/2)WD^(−1/2): normalized adjacency matrix of G
    d_negative_half = sparse_d.power(-0.5)
    sparse_w_tilde = d_negative_half @ sparse_w @ d_negative_half
    #L˜ = D − W˜ : Laplacian of G
    sparse_l_tilde = sparse_d - sparse_w_tilde
    #creating I indentity matrix
    identity = [1] * len(list(g.nodes))
    identity_index = list(range(0, node_len))
    sparse_i = csc_matrix((identity, (identity_index, identity_index)),shape=(node_len,node_len))
    #(I + alpha*L~)^-1
    almost_reg_laplacian = sparse_i + (alpha * sparse_l_tilde)
    reg_laplacian = inv(almost_reg_laplacian)
    
    return reg_laplacian

def regularized_laplacian(go_id):
    g = human_ppi_parser()
    reg_laplacian = rl(g,1)

    #s(u) = (I + alpha*L~)^-1 * y(u) [y -> 1/nodes_visited:0 positive:otherwise]
    #s(u) = reg_laplacian * y