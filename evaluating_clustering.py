import sys
import os

import pandas as pd
import csv

import spaligner_parser


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

def evaluate_clustering(clustering_tsv, spaligner_clustering_tsv):
    clusters = tsv_to_sets(clustering_tsv)
    spaligner_clusters = tsv_to_sets(spaligner_clustering_tsv)

    J = jaccard_similarity(clusters, spaligner_clusters)
    print('Jaccard similarity: %.3f' % J)


def main():
    outdir = os.path.join(sys.argv[3])
    spaligner_clustering_tsv = spaligner_parser.spaligner_to_clustering_tsv(sys.argv[2], os.path.join(outdir, 'spaligner_clustering.tsv'))
    evaluate_clustering(sys.argv[1], spaligner_clustering_tsv)


if __name__ == '__main__':
    main()