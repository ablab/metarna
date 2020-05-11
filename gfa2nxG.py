# gfa2nxG.py assembly_graph_with_scaffolds.gfa

import sys

import scipy as sp

import networkx as nx

def line_to_node(line):
    fields = line.strip().split()
    name = fields[1]
    attr = {'len': len(fields[2])}
    if 'KC:i:' in line:
        kmer_count = int(fields[3][5:])
        attr['cov'] = kmer_count
    return name, attr

def line_to_edge(line):
    fields = line.strip().split()
    u = fields[1]
    v = fields[3]
    attr = {'FromOrient': fields[2], 'ToOrient': fields[4], 'cigar': fields[5]}
    return u, v, attr


def gfa_to_G(gfa):
    G = nx.MultiDiGraph()
    with open(gfa, 'r') as fin:
        for line in fin:
            record_type = line[0]
            if record_type in ['#', 'H', 'C', 'P']:
                continue
            elif record_type == 'S':
                name, attr = line_to_node(line)
                G.add_node(name, **attr)
            elif record_type == 'L':
                u, v, attr = line_to_edge(line)
                G.add_edge(u, v, **attr)
    return G

def get_A(G):
    A = nx.adjacency_matrix(G)
    # print(A.todense())
    return A

def get_X(G):
    X = []
    for node in G.nodes:
        X.append(list(G.nodes[node].values()))
    # print(X)
    return X

# for node in G:
    # for nbrs in G[node]:
    #     for edge in G[node][nbrs]:
    #         print(list(G[node][nbrs][edge].values()))

gfa = sys.argv[1]

G = gfa_to_G(gfa)
A = get_A(G)
X = get_X(G)