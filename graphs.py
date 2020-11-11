import networkx as nx

from collections import defaultdict

import itertools

import numpy as np

from scipy.stats.mstats import gmean
from scipy.stats import hmean

import time

import sys
import os

import matplotlib.pyplot as plt

from spaligner_parser import spaligner_to_df_not_ss


def filter_G_by_degree(G, filtered_degree=2):
    removed = [node for node, degree in dict(G.degree()).items() if degree < filtered_degree]
    G.remove_nodes_from(removed)
    return G

def filter_G_by_weight(G, weight, treshold):
    weights = nx.get_edge_attributes(G, weight)
    w_list = list(weights.values())
    k = max(0, int(len(w_list) * treshold) - 1)
    min_weight = np.partition(w_list, k)[k]
    G.remove_edges_from([e for e, w in weights.items() if w < min_weight])
    G.name += '_min_weight_{:.3f}'.format(min_weight)
    write_G_statistics(G)

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

def get_weight_attr(G, u, v, num_long_reads=0):
    cov = nx.get_node_attributes(G, 'cov')
    cov_diff = 1.0 / (abs(cov[u] - cov[v]) + sys.float_info.epsilon)
    weight_attr = {'cov_diff': cov_diff,
                   'long_reads': num_long_reads,
                   'geometric_mean': gmean([cov_diff, num_long_reads]),
                   'harmonic_mean': hmean([cov_diff, num_long_reads + sys.float_info.epsilon])}
    return weight_attr

def write_G_statistics(G):
    print('{} graph statistics: Nodes: {}; Edges: {}; Components: {}.'.
          format(G.name, G.number_of_nodes(),
                 G.number_of_edges(),
                 nx.number_connected_components(G)))

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
    # edges = set()
    weight_attr = {}
    num_long_reads = defaultdict(int)
    start = time.time()
    tsv_df = spaligner_to_df_not_ss(spaligner_tsv, G)
    for path_str in tsv_df['path of the alignment']:
        path = path_str.replace(';', ',').split(',')
        for u, v in itertools.combinations(path, 2):
            # edges.add((u, v))
            num_long_reads[(u, v)] += 1
            weight_attr[(u, v)] = get_weight_attr(G, u, v, num_long_reads[(u, v)])
    end = time.time()
    # print('Elapsed time on long reads graph construction: {}'.format((end - start) * 1.0 / 60 / 60))
    return weight_attr

def G_to_friendships_graph(G, spaligner_long_reads_tsv):
    # fG = G.to_undirected()
    fG = G.copy()
    fG.name = 'friendships'
    # fG.add_edges_from(gfa_parser.get_friendships(G))
    weight_attr = get_friendships_from_long_reads(spaligner_long_reads_tsv, fG)
    fG.add_edges_from((edge[0], edge[1], w_dict) for edge, w_dict in weight_attr.items())
    write_G_statistics(fG)
    return fG

def truncate_values(w_dict, keys):
    truncated_dict = {}
    for key, value in w_dict.items():
        if key in keys:
            truncated_dict[key] = str(value)[:5]
    return truncated_dict

def plot_graph_components(G, weight, outdir, n=4):
    options = {'with_labels': True,
               'pos': nx.spring_layout(G),
               'font_size': 5}
    largest_components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[:n]
    for i, component in enumerate(largest_components):
        nx.draw(component, **options)
        nx.draw_networkx_edge_labels(component, options['pos'], font_size=options['font_size'],
                                     edge_labels=truncate_values(nx.get_edge_attributes(G, weight), component.edges))
        plt.savefig(os.path.join(outdir, '{}.{}.png'.format(G.name, i)))
        plt.clf()