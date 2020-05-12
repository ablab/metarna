# gfa2nxG.py assembly_graph_with_scaffolds.gfa alignment.tsv

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
    e_key = (fields[2], fields[4])
    attr = {'cigar': fields[5]}
    return u, v, e_key, attr

def set_edge_features(G):
    from Bio.Seq import reverse_complement

    for n, nbrsdict in G.adjacency():
        for nbr, edict in nbrsdict.items():
            for e_key, eattr in edict.items():
                seq_n = G.nodes[n]['seq']
                if e_key[0] == '-':
                    seq_n = reverse_complement(seq_n)
                seq_nbr = G.nodes[nbr]['seq']
                if e_key[1] == '-':
                    seq_nbr = reverse_complement(seq_nbr)
                overlap = int(G.edges[n, nbr, e_key]['cigar'][:-1])
                G.edges[n, nbr, e_key]['seq'] = seq_n + seq_nbr[overlap:]
                G.edges[n, nbr, e_key]['len'] = len(G.edges[n, nbr, e_key]['seq'])
                G.edges[n, nbr, e_key]['cov'] = G.nodes[n]['KC'] + G.nodes[nbr]['KC'] - 1
                for nucl in ['A', 'C', 'G', 'T']:
                    G.edges[n, nbr, e_key][nucl] = G.edges[n, nbr, e_key]['seq'].count(nucl)
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
                u, v, e_key, attr = line_to_edge(line)
                G.add_edge(u, v, key=e_key, **attr)
    return G

def get_A(G):
    A = nx.adjacency_matrix(G)
    # print(A.todense())
    return A

def get_X(G):
    X = []
    for n, nbrsdict in G.adjacency():
        for nbr, edict in nbrsdict.items():
            for e_key, eattr in edict.items():
                X.append([eattr[key] for key in ['len', 'cov', 'A', 'C', 'G', 'T']])
    # print(X)
    return X

def path_to_edges(path):
    # path is something like 123+;288-,128+
    edges = []
    alignments = path.split(';')
    for alignment in alignments:
        records = alignment.split(',')
        for i in range(len(records) - 1):
            u = records[i][:-1]
            v = records[i + 1][:-1]
            e_key = (records[i][-1], records[i + 1][-1])
            edges.append((u, v, e_key))
    return edges

def set_edge_labels(G, tsv):
    import pandas as pd

    tsv_df = pd.read_csv(tsv, sep="\t", names=['sequence name',
                                               'start position of alignment on sequence',
                                               'end position of  alignment on sequence',
                                               'start position of alignment on the first edge of the Path',
                                               'end position of alignment on the last edge of the Path',
                                               'sequence length',
                                               'path of the alignment',
                                               'lengths of the alignment on each edge of the Path respectively',
                                               'sequence of alignment Path'])

    # Split path column into multiple rows
    new_df = pd.DataFrame(tsv_df['path of the alignment'].apply(path_to_edges).tolist(), index=tsv_df['sequence name']).stack()
    new_df = new_df.reset_index([0, 'sequence name'])
    new_df.columns = ['labels', 'edge']

    # Generate list of sequence names for each edge
    grouped_df = new_df.groupby('edge')['labels'].apply(list).reset_index()

    nx.set_edge_attributes(G, grouped_df.set_index('edge')['labels'].to_dict(), name='labels')
    # print(nx.get_edge_attributes(G,'labels'))
    return G

def set_node_labels(G, tsv):
    import pandas as pd

    tsv_df = pd.read_csv(tsv, sep="\t", names=['sequence name',
                                               'start position of alignment on sequence',
                                               'end position of  alignment on sequence',
                                               'start position of alignment on the first edge of the Path',
                                               'end position of alignment on the last edge of the Path',
                                               'sequence length',
                                               'path of the alignment',
                                               'lengths of the alignment on each edge of the Path respectively',
                                               'sequence of alignment Path'])

    # Split path column into multiple rows
    new_df = pd.DataFrame(tsv_df['path of the alignment'].str.replace(';', ',').str.split(',').tolist(),
                          index=tsv_df['sequence name']).stack()
    new_df = new_df.reset_index([0, 'sequence name'])
    new_df.columns = ['sequence name', 'node']

    # Generate list of sequence names for each node with orientation
    grouped_df = new_df.groupby('node')['sequence name'].apply(list).reset_index()

    grouped_dict = grouped_df.set_index('node')['sequence name'].to_dict()
    labels_dict = {}
    for node in G.nodes:
        forward_label = grouped_dict[node + '+'] if node + '+' in grouped_dict else []
        reverse_label = grouped_dict[node + '-'] if node + '-' in grouped_dict else []
        labels_dict[node] = (forward_label, reverse_label)
    nx.set_node_attributes(G, labels_dict, name='label')
    # print(nx.get_node_attributes(G,'label'))
    return G

# SPAdes output
gfa = sys.argv[1]

# SPAligner output
tsv = sys.argv[2]

# Get graph from gfa file
G = gfa_to_G(gfa)

# Get something like Adjacency matrix
# The weights are summed for MultiDiGraph parallel edges
A = get_A(G)

# Get feature matrix
G = set_edge_features(G)
X = get_X(G)

# Get labels for edges
G = set_edge_labels(G, tsv)
