import networkx as nx
from scipy.sparse import csc_matrix, identity, eye
from scipy.sparse.linalg import inv
import scipy.sparse as sp
import pandas as pd
import os

from parser_module import Annotation_Parser, ROOT_TERMS
import numpy as np

node_index = {}

def human_ppi_parser():
    df = pd.read_csv("human_interactome_weighted.txt", sep="\t")
    g = nx.Graph()
    index = 0

    for ind in df.index:
        data = [df['#tail'][ind], df['head'][ind]]
        if data[0] not in g.nodes():
            g.add_node(data[0])
            node_index[data[0]] = index
            index += 1
        if data[1] not in g.nodes():
            g.add_node(data[1])
            node_index[data[1]] = index
            index += 1
        
        g.add_edge(data[1], data[0])

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

        row.append(node_index[edge[1]])
        coln.append(node_index[edge[0]])
        data.append(1)

    sparse_w = csc_matrix((data, (row, coln)),shape=(node_len,node_len))

    #D: diagonal matrix where Duu = Sum[v] w_uv , for every node u in G
    flattened_sum_w = np.asarray(sparse_w.sum(axis=0)).flatten()
    flattened_sum_w[np.isinf(flattened_sum_w)] = 0
    sparse_d = sp.diags(flattened_sum_w)

    #W˜ = D^(−1/2)WD^(−1/2): normalized adjacency matrix of G
    d_negative_half = sparse_d.power(-0.5)
    sparse_w_tilde = (d_negative_half @ sparse_w) @ d_negative_half

    flattened_sum_w_tilde = np.asarray(sparse_w_tilde.sum(axis=0)).flatten()
    flattened_sum_w_tilde[np.isinf(flattened_sum_w_tilde)] = 0
    sparse_d = sp.diags(flattened_sum_w_tilde)

    #L˜ = D − W˜ : Laplacian of G
    sparse_l_tilde = sparse_d - sparse_w_tilde
    #creating I indentity matrix
    sparse_i = identity(node_len)
    #(I + alpha*L~)^-1
    almost_reg_laplacian = sparse_i + (alpha * sparse_l_tilde)
    reg_laplacian = inv(almost_reg_laplacian)

    return reg_laplacian

#will run for each go_term, printing out results in their own file within rl_output folder
def all_regularized_laplacian(alpha):

    g = human_ppi_parser()
    reg_laplacian = rl(g,alpha)

    #default threshold is 1%
    #make sure these files are in the expected directory
    data = Annotation_Parser('go.obo', 'goa_human.gaf', 'human_interactome_weighted.txt', threshold=10)
    #data.stacked_plot()
    node_len = len(list(g.nodes))
    for term in ["GO:0048523","GO:0019538","GO:0043167","GO:0005515","GO:0005576","GO:0005829"]:
        pos,neg = data.examples(term)
        val = []
        row = []
        coln = [0] * node_len
        for x in list(g.nodes()):
            val.append(pos[data.gene_index(x)])
            row.append(node_index[x])

        outfile = open(f"rl outputs/{term[3:]}.txt", 'w')
        outfile.write("Node  Name  Score  Value  Example\n")
        y = csc_matrix((val, (row, coln)),shape=(node_len,1))
        s = reg_laplacian @ y
        for i in node_index:
            example = "Positive Example" if (int(y[node_index[i]].toarray()[0][0]) == 1) else "Non Example"
            outfile.write(f"Node:  {i}  Score:  {s[node_index[i]].toarray()[0][0]}  {example}\n")
        outfile.close()
        
def murali_group(alpha):
    """
    Normalize W by multiplying D^(-1/2) * W * D^(-1/2)
    This is used for GeneMANIA
    *W*: weighted network as a scipy sparse matrix in csr format
    """
    g = human_ppi_parser()
    #W : adjacency matrix of G
    node_len = len(list(g.nodes))
    row = []
    coln = []
    data = []
    for edge in g.edges():
        row.append(node_index[edge[0]])
        coln.append(node_index[edge[1]])
        data.append(1)

        row.append(node_index[edge[1]])
        coln.append(node_index[edge[0]])
        data.append(1)
    
    W = csc_matrix((data, (row, coln)),shape=(node_len,node_len))

    deg = np.asarray(W.sum(axis=0)).flatten()
    deg = np.divide(1., np.sqrt(deg))
    deg[np.isinf(deg)] = 0
    D = sp.diags(deg)
    # normalize W by multiplying D^(-1/2) * W * D^(-1/2)
    P = D.dot(W.dot(D))

    deg = np.asarray(P.sum(axis=0)).flatten()
    deg[np.isinf(deg)] = 0
    D = sp.diags(deg)
    L = D - P

    M = eye(L.shape[0]) + (alpha*L)     
    # now take the inverse
    m_inv = inv(M)

    #default threshold is 1%
    #make sure these files are in the expected directory
    data = Annotation_Parser('go.obo', 'goa_human.gaf', 'human_interactome_weighted.txt', threshold=10)
    #data.stacked_plot()
        
    for term in ["GO:0048523","GO:0019538","GO:0043167","GO:0005515","GO:0005576","GO:0005829"]:
        pos,neg = data.examples(term)
        val = []
        row = []
        coln = [0] * node_len
        for x in list(g.nodes()):
            val.append(pos[data.gene_index(x)])
            row.append(node_index[x])

        outfile = open(f"rl outputs/muraligroup_{term[3:]}.txt", 'w')
        outfile.write("Node  Name  Score  Value  Example\n")
        y = csc_matrix((val, (row, coln)),shape=(node_len,1))
        s = m_inv.dot(y)
        for i in node_index:
            example = "Positive Example" if (int(y[node_index[i]].toarray()[0][0]) == 1) else "Non Example"
            outfile.write(f"Node:  {i}  Score:  {s[node_index[i]].toarray()[0][0]}  {example}\n")
        outfile.close()

all_regularized_laplacian(1)
murali_group(1)
 
# iterate over files in directory
for filename in os.listdir("rl outputs"):
    outfile = os.path.join("rl outputs", filename)
    # checking if it is a file
    if os.path.isfile(outfile):
        dataframe = pd.read_csv(outfile, sep="  ")
        dataframe = dataframe.sort_values(by=['Value'], ascending=False)
        dataframe.to_csv(outfile, sep='\t', index=False)
        out = open(outfile, 'r')
        final = []
        lines = out.readlines()

        skip = True
        for line in lines:
            if not skip:
                final.append(line.replace('\t', '  '))
            else:
                skip = False
        out.close()
        out = open(outfile, 'w')
        for line in final:
            out.write(line)
        out.close()