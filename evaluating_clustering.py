import sys
import os

import pandas as pd
import csv

import spaligner_parser
from gfa_parser import gfa_to_G

def tsv_to_sets(tsv):
    clusters = set()
    with open(tsv, 'r') as fin:
        for line in fin:
            clusters.add(frozenset(line.split()))
    return clusters

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

def F1_best_match(r_cluster, ground_truth_set):
    F1_best_match = 0
    for gt_cluster in ground_truth_set:
        F1_curr = F1_for_two_clusters(r_cluster, gt_cluster)
        if F1_best_match < F1_curr:
            F1_best_match = F1_curr
    return F1_best_match

def F1_for_clustering(reconstructed_set, ground_truth_set):
    F1 = 0
    for r_cluster in reconstructed_set:
        F1 += F1_best_match(r_cluster, ground_truth_set)
    F1 /= len(reconstructed_set)
    return F1

def evaluate_clustering(reconstructed_clustering_tsv, ground_truth_clustering_tsv):
    reconstructed_clusters = tsv_to_sets(reconstructed_clustering_tsv)
    ground_truth_clusters = tsv_to_sets(ground_truth_clustering_tsv)

    J = jaccard_similarity(reconstructed_clusters, ground_truth_clusters)
    print('Jaccard similarity: %.3f' % J)

    F1 = F1_for_clustering(reconstructed_clusters, ground_truth_clusters)
    print('F1 score: %.3f' % F1)


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
    evaluate_clustering(clustering_tsv, spaligner_clustering_tsv)


if __name__ == '__main__':
    main()