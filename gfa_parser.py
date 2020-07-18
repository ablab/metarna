# gfa_parser.py assembly_graph_with_scaffolds.gfa alignment.tsv outdir

import sys, os

import time

import scipy as sp
import pandas as pd

import networkx as nx

from Bio.Seq import reverse_complement

from spaligner_parser import spaligner_to_df_not_ss


def get_one_type_gfa(gfa, type, outdir):
    one_type_gfa = os.path.join(outdir, '{}.gfa'.format(type))
    os.system('grep \'^{}\' {} > {}'.format(type, gfa, one_type_gfa))
    return one_type_gfa

# 'RecordType', 'PathName', 'SegmentNames', 'Overlaps'
def one_type_gfa_to_df(one_type_gfa):
    p_df = pd.read_csv(one_type_gfa,
                       sep="\t", header=None, usecols=[1, 2],
                       names=['PathName', 'SegmentNames'])
    return p_df

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

def get_X(nodes, out_tsv):
    X = []
    features = ['len', 'cov', 'A', 'C', 'G', 'T']
    with open(out_tsv, 'w') as fout:
        fout.write('node ' + ' '.join(features) + '\n')
        for node in nodes:
            X.append([nodes[node][key] for key in features])
            fout.write(node + ' ' + str(X[-1][0]) + ' ' + ' '.join(["%.2f" % e for e in X[-1][1:]]) + '\n')
    return X

# Path existence between nodes in gfa graph means edge in friendship graph
# Nodes connected in friendship graph more likely to belong one gene (like social community)
def get_friendships(G):
    start = time.time()
    friendships = [(u, v) for u in G.nodes for v in G.nodes
                   if nx.algorithms.shortest_paths.generic.has_path(G, u, v)]
    end = time.time()
    print('Elapsed time on friendship graph construction: ', (end - start) * 1.0 / 60 / 60)
    return friendships

def get_friendships_from_long_reads(spaligner_tsv, G):
    friendships = []
    start = time.time()
    tsv_df = spaligner_to_df_not_ss(spaligner_tsv, G)
    for path_str in tsv_df['path of the alignment']:
        path = path_str.replace(';', ',').split(',')
        friendships += [(u, v) for i, u in enumerate(path)
                        for j, v in enumerate(path) if i < j]
    end = time.time()
    print('Elapsed time on long reads graph construction: ', (end - start) * 1.0 / 60 / 60)
    return friendships

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
    X = get_X(G.nodes, features_tsv)


if __name__ == '__main__':
    main()