# gfa2nxG.py assembly_graph_with_scaffolds.gfa

import sys

import scipy as sp

import networkx as nx

def line_to_node(line):
    fields = line.strip().split()
    name = fields[1]
    attr = {'seq': fields[2]}
    if 'KC:i:' in line:
        kmer_count = int(fields[3][5:])
        attr['KC'] = kmer_count
    return name, attr

def line_to_edge(line):
    fields = line.strip().split()
    u = fields[1]
    v = fields[3]
    attr = {'FromOrient': fields[2], 'ToOrient': fields[4], 'cigar': fields[5]}
    return u, v, attr

def set_edge_features(G):
    from Bio.Seq import reverse_complement

    for n, nbrsdict in G.adjacency():
        for nbr, edict in nbrsdict.items():
            for e, eattr in edict.items():
                seq_n = G.nodes[n]['seq']
                if G.edges[n, nbr, e]['FromOrient'] == '-':
                    seq_n = reverse_complement(seq_n)
                seq_nbr = G.nodes[nbr]['seq']
                if G.edges[n, nbr, e]['ToOrient'] == '-':
                    seq_nbr = reverse_complement(seq_nbr)
                overlap = int(G.edges[n, nbr, e]['cigar'][:-1])
                G.edges[n, nbr, e]['seq'] = seq_n + seq_nbr[overlap:]
                G.edges[n, nbr, e]['len'] = len(G.edges[n, nbr, e]['seq'])
                G.edges[n, nbr, e]['cov'] = G.nodes[n]['KC'] + G.nodes[nbr]['KC'] - 1
                for nucl in ['A', 'C', 'G', 'T']:
                    G.edges[n, nbr, e][nucl] = G.edges[n, nbr, e]['seq'].count(nucl)
    return G

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
    for n, nbrsdict in G.adjacency():
        for nbr, edict in nbrsdict.items():
            for e, eattr in edict.items():
                X.append([eattr[key] for key in ['len', 'cov', 'A', 'C', 'G', 'T']])
    # print(X)
    return X

gfa = sys.argv[1]

# Get graph from gfa file
G = gfa_to_G(gfa)

# Get something like Adjacency matrix
# The weights are summed for MultiDiGraph parallel edges
A = get_A(G)

# Get feature matrix
G = set_edge_features(G)
X = get_X(G)

