# gfa2nxG.py assembly_graph_with_scaffolds.gfa alignment.tsv outdir

import sys, os

import scipy as sp
import pandas as pd

import networkx as nx

from Bio.Seq import reverse_complement


def line_to_node(line):
    fields = line.strip().split()
    name = fields[1]
    attr = {'seq': fields[2]}
    if 'KC:i:' in line:
        kmer_count = int(fields[3][5:])
        attr['KC'] = kmer_count
    return name, attr

# L       934049  -       36137   +       49M
def line_to_edge(line):
    fields = line.strip().split()
    # Node name plus node orientation
    u = fields[1] + fields[2]
    v = fields[3] + fields[4]
    attr = {'cigar': fields[5]}
    return u, v, attr

def line_to_rc_edge(line):
    rc_dict = {'+': '-', '-': '+'}
    fields = line.strip().split()
    # Node name plus node orientation
    u = fields[3] + rc_dict[fields[4]]
    v = fields[1] + rc_dict[fields[2]]
    attr = {'cigar': fields[5]}
    return u, v, attr

def gfa_to_G(gfa):
    G = nx.DiGraph()
    with open(gfa, 'r') as fin:
        for line in fin:
            record_type = line[0]
            if record_type in ['#', 'H', 'C', 'P']:
                continue
            elif record_type == 'S':
                name, attr = line_to_node(line)
                G.add_node(name + '+',
                           seq=attr['seq'],
                           cov=attr['KC'] * 1.0 / len(attr['seq']),
                           len=len(attr['seq']),
                           A=attr['seq'].count('A') * 1.0 / len(attr['seq']),
                           C=attr['seq'].count('C') * 1.0 / len(attr['seq']),
                           G=attr['seq'].count('G') * 1.0 / len(attr['seq']),
                           T=attr['seq'].count('T') * 1.0 / len(attr['seq']))
                G.add_node(name + '-',
                           seq=reverse_complement(attr['seq']),
                           cov=attr['KC'] * 1.0 / len(attr['seq']),
                           len=len(attr['seq']),
                           A=attr['seq'].count('T') * 1.0 / len(attr['seq']),
                           C=attr['seq'].count('G') * 1.0 / len(attr['seq']),
                           G=attr['seq'].count('C') * 1.0 / len(attr['seq']),
                           T=attr['seq'].count('A') * 1.0 / len(attr['seq']))
            elif record_type == 'L':
                u, v, attr = line_to_edge(line)
                G.add_edge(u, v, **attr)

                u, v, attr = line_to_rc_edge(line)
                G.add_edge(u, v, **attr)
    return G

def get_A(G):
    A = nx.adjacency_matrix(G)
    print(A.todense())
    return A

def get_X(G, out_tsv):
    X = []
    features = ['len', 'cov', 'A', 'C', 'G', 'T']
    with open(out_tsv, 'w') as fout:
        fout.write('node ' + ' '.join(features) + '\n')
        for node in G.nodes:
            X.append([G.nodes[node][key] for key in features])
            fout.write(node + ' ' + str(X[-1][0]) + ' ' + ' '.join(["%.2f" % e for e in X[-1][1:]]) + '\n')
    return X

def spaligner_to_df(tsv):
    tsv_df = pd.read_csv(tsv, sep="\t", names=['sequence name',
                                               'start position of alignment on sequence',
                                               'end position of  alignment on sequence',
                                               'start position of alignment on the first edge of the Path',
                                               'end position of alignment on the last edge of the Path',
                                               'sequence length',
                                               'path of the alignment',
                                               'lengths of the alignment on each edge of the Path respectively',
                                               'sequence of alignment Path'])
    return tsv_df

def set_node_labels(G, tsv, node_to_db_tsv):
    tsv_df = spaligner_to_df(tsv)

    # Split path column into multiple rows
    new_df = pd.DataFrame(tsv_df['path of the alignment'].str.replace(';', ',').str.split(',').tolist(),
                          index=tsv_df['sequence name']).stack()
    new_df = new_df.reset_index([0, 'sequence name'])
    new_df.columns = ['sequence name', 'node']

    # Generate set of sequence names for each node with orientation
    grouped_df = new_df.groupby('node')['sequence name'].apply(set).reset_index()

    grouped_dict = grouped_df.set_index('node')['sequence name'].to_dict()
    nx.set_node_attributes(G, grouped_dict, name='label')

    # print(nx.get_node_attributes(G,'label'))

    with open(node_to_db_tsv, 'w') as fin:
        for node, transcripts in grouped_dict.items():
            fin.write(node + '\t' + ' '.join(transcripts) + '\n')

    return G

def main():
    # SPAdes output
    gfa = sys.argv[1]

    # SPAligner output
    tsv = sys.argv[2]

    outdir = sys.argv[3]

    # Get graph from gfa file
    G = gfa_to_G(gfa)

    # Get Adjacency matrix
    A = get_A(G)

    # Get feature matrix
    features_tsv = os.path.join(outdir, 'features.tsv')
    X = get_X(G, features_tsv)

    # Set labels for nodes
    node_to_db_tsv = os.path.join(outdir, 'node_to_db.tsv')
    G = set_node_labels(G, tsv, node_to_db_tsv)


if __name__ == '__main__':
    main()