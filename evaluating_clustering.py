import sys
import os

import pandas as pd
import csv

import matplotlib.pyplot as plt
plt.switch_backend('agg')

import networkx as nx

import spaligner_parser
from gfa_parser import gfa_to_G

import graphs
import clustering


def tsv_to_sets(tsv, min_component_size=3):
    clusters = set()
    with open(tsv, 'r') as fin:
        for line in fin:
            path = frozenset(line.strip().split(','))
            if len(path) < min_component_size:
                continue
            clusters.add(path)
    print(tsv + ': {} clusters'.format(len(clusters)))
    return clusters


def write_clustering(clustering, tsv, min_clusters_size=2):
    with open(tsv, 'w') as outfile:
        for cluster in clustering:
            if len(cluster) < min_clusters_size:
                continue
            outfile.write(','.join([str(x) for x in cluster]) + '\n')


def jaccard_similarity(set1, set2):
    up = len(set1.intersection(set2))
    down = len(set1.union(set2))
    # print('Intersection:')
    # for c in set1.intersection(set2):
    #     print(c)
    # print('Ground truth - clustering: ')
    # for c in set2 - set1:
    #     print(c)
    print('Exact reconstruction: {}'.format(up))
    print('Clusters in total: {}'.format(down))
    return up / down

def F1_for_two_clusters(reconstructed_cluster, ground_truth_cluster):
    precision = \
        len(ground_truth_cluster.intersection(reconstructed_cluster)) / \
        len(reconstructed_cluster)
    recall = \
        len(ground_truth_cluster.intersection(reconstructed_cluster)) / \
        len(ground_truth_cluster)
    if precision + recall != 0:
        F1 = 2 * precision * recall / (precision + recall)
    else:
        F1 = 0
    return F1

def F1_best_match(r_cluster, ground_truth_set, fout):
    F1_best_match = 0
    cluster_best_match = None
    for gt_cluster in ground_truth_set:
        F1_curr = F1_for_two_clusters(r_cluster, gt_cluster)
        if F1_best_match <= F1_curr:
            F1_best_match = F1_curr
            cluster_best_match = gt_cluster
    if F1_best_match != 1:
        fout.write(' '.join(sorted(r_cluster)) + '\n' + ' '.join(sorted(cluster_best_match)) + '\n\n')
    return F1_best_match

def F1_for_clustering(reconstructed_set, ground_truth_set, outdir):
    F1 = 0
    not_reconstructed_txt = os.path.join(outdir, 'not_reconstructed.debug')
    with open(not_reconstructed_txt, 'w') as fout:
        for r_cluster in reconstructed_set:
            F1 += F1_best_match(r_cluster, ground_truth_set, fout)
        F1 /= len(reconstructed_set)
    return F1

def exact_recall(reconstructed_set, ground_truth_set):
    up = len(reconstructed_set.intersection(ground_truth_set))
    down = len(ground_truth_set)
    return up / down

def evaluate_clustering(reconstructed_clustering_tsv, ground_truth_clustering_tsv, outdir):
    short_report_txt = os.path.join(outdir, 'short_report.txt')

    reconstructed_clusters = tsv_to_sets(reconstructed_clustering_tsv)
    ground_truth_clusters = tsv_to_sets(ground_truth_clustering_tsv)

    J = jaccard_similarity(reconstructed_clusters, ground_truth_clusters)
    print('Jaccard similarity: %.3f' % J)

    recall = exact_recall(reconstructed_clusters, ground_truth_clusters)
    print('Recall: %.3f' % recall)

    F1 = F1_for_clustering(reconstructed_clusters, ground_truth_clusters, outdir)
    print('F1 score: %.3f' % F1)

    with open(short_report_txt, 'w') as fout:
        fout.write('Jaccard similarity: %.3f\n' % J)
        fout.write('Recall: %.3f\n' % recall)
        fout.write('F1 score: %.3f\n' % F1)

def get_node_colors(G, c_dict):
    clusters = []
    for node in G.nodes:
        clusters.append(c_dict[node])
    size = len(set(clusters))
    node_colors = [str(cluster * 1.0 / size) for cluster in clusters]
    return node_colors

def plot_components_clusters(G, c_list, weight, outdir, n=4):
    c_dict = clustering.clusters_list_to_dict(c_list)
    pos = nx.spring_layout(G)
    largest_components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[:n]
    for i, component in enumerate(largest_components):
        colors = get_node_colors(component, c_dict)
        edge_labels = graphs.truncate_values(nx.get_edge_attributes(G, weight), component.edges)
        nx.draw_networkx_nodes(component, pos=pos, node_color=colors)
        nx.draw_networkx_labels(component, pos, font_size=5)
        nx.draw_networkx_edges(component, pos, alpha=0.5)
        nx.draw_networkx_edge_labels(component, pos=pos, font_size=5, edge_labels=edge_labels)
        plt.savefig(os.path.join(outdir, '{}.{}.png'.format(G.name, i)))
        plt.clf()

def plot_graph_clusters(G, c_list, outdir):
    size = float(len(c_list))
    pos = nx.spring_layout(G)
    for i, com in enumerate(c_list):
        nx.draw_networkx_nodes(G, pos, com, node_size=20, node_color=str(i / size))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.savefig(os.path.join(outdir, '{}.png'.format(G.name)))
    plt.clf()


def main():
    clustering_tsv = sys.argv[1]
    spaligner_tsv = sys.argv[2]
    gfa = sys.argv[3]
    k = int(sys.argv[4])
    outdir = os.path.join(sys.argv[5])

    spaligner_clustering_tsv = \
        spaligner_parser.spaligner_to_clustering_tsv(spaligner_tsv,
                                                     os.path.join(outdir, 'spaligner_clustering.tsv'),
                                                     gfa_to_G(gfa, k))
    evaluate_clustering(clustering_tsv, spaligner_clustering_tsv, outdir)


if __name__ == '__main__':
    main()