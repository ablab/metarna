# nxG2clusters.py assembly_graph_with_scaffolds.gfa outdir

import sys
import os

import time

import argparse

import networkx as nx

import matplotlib.pyplot as plt
plt.switch_backend('agg')

import pandas as pd

from sklearn.preprocessing import StandardScaler

from persona.persona import CreatePersonaGraph
from persona.directed_persona import CreateDirectedPersonaGraph
from persona.persona import PersonaOverlappingClustering
from persona.flags import _CLUSTERING_FN
from persona.splitter import do_embedding

from gfa_parser import gfa_to_G
import graphs
import spaligner_parser
import visualising_embedding
import evaluating_clustering


def parse_args():
    parser = argparse.ArgumentParser(description='Clustering on graphs',
                                     usage='{} --gfa assembly_graph_with_scaffolds.gfa '
                                           '--ground_truth transcripts_alignment.tsv '
                                           '--friendships_reads reads_alignment.tsv '
                                           '-k 49 --outdir clustering_out'.format(sys.argv[0]))
    parser.add_argument('--clustering', '-c', dest='c_name', default='best_partition',
                        required=False, type=str,
                        help='Choose the algorithm for local and global clustering',
                        choices=['label_prop', 'modularity', 'connected_components',
                                 'weakly_connected_components', 'best_partition'])
    parser.add_argument('--weight', '-w', dest='w_name', default='reads_and_db',
                        required=False, type=str,
                        help='Choose the weight for clustering',
                        choices=['cov_diff', 'reads_and_db', 'geometric_mean', 'harmonic_mean'])
    parser.add_argument('--gfa', '-g', required=True, help='Assembly graph')
    parser.add_argument('--ground_truth', dest='spaligner_ground_truth_tsv', required=True,
                        help='It can be transcripts aligned to assembly graph using SPAligner [tsv]',)
    parser.add_argument('--friendships_reads', dest='spaligner_long_reads_tsv', required=False,
                        help='Long reads aligned to assembly graph '
                             '(or any other confirmation of belonging to one transcript) [tsv]')
    parser.add_argument('--friendships_db', dest='spaligner_db_tsv', required=False,
                        help='Reference transcripts aligned to assembly graph '
                             '(or any other confirmation of belonging to one transcript) [tsv]')
    parser.add_argument('-k', type=int, required=True,
                        help='k-mer value used in assembly graph construction')
    parser.add_argument('--outdir', '-o', required=True)
    parser.add_argument('--filter', default=None, type=float,
                        help='Filter this percent of edges based on their weights')

    args = parser.parse_args()

    return args


def remove_regular_model(in_path, out_path):
    fout = open(out_path, 'w')

    with open(in_path, 'r') as fin:
        for line in fin:
            node = line.split()[0]
            if '+_' in node or '-_' in node:
                fout.write(line)

    fout.close()
    return out_path


def get_tst_G(G):
    # path1 930004-,278546-,36185+,278990+,283130+,352975-,37703+
    # path2 930004-,239212-,36185+,365256-,283130+,352975-,37703+
    nodes_tst = ['36185+', '37703+', '239212-', '278546-', '278990+',
                 '283130+', '352975-', '365256-', '930004-', '2326645-']
    G_tst = G.subgraph(nodes_tst).copy()
    return G_tst


def get_total_emb(p_emb_tsv, features_tsv, persona_to_node_tsv):
    # concatenate structural features (persona graph embedding)
    # and node features (len, cov, A, C, G, T)
    p_emb = pd.read_csv(p_emb_tsv, sep=' ', header=None, index_col=0)

    features_df = pd.read_csv(features_tsv, sep=' ',
                           header=None, index_col=0, skiprows=1,
                           names=range(p_emb.shape[1], p_emb.shape[1] + 7))
    # It will be helpful to convert each feature into z-scores
    # (number of standard deviations from the mean) for comparability
    scaled_features = StandardScaler().fit_transform(features_df.values)
    scaled_features_df = pd.DataFrame(scaled_features, index=features_df.index, columns=features_df.columns)

    persona_to_node = pd.read_csv(persona_to_node_tsv, sep=' ',
                                  header=None, index_col=0,
                                  names=['initial_node'])

    tot_emb_df = pd.concat([p_emb, persona_to_node], axis=1).join(scaled_features_df, on='initial_node')
    tot_emb_df = tot_emb_df.drop(columns=['initial_node'])

    return tot_emb_df


def main():
    args = parse_args()

    local_clustering_fn = _CLUSTERING_FN[args.c_name]
    global_clustering_fn = _CLUSTERING_FN[args.c_name]

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    G = gfa_to_G(args.gfa, args.k)

    # G = get_tst_G(G)

    # G = graphs.filter_G_by_degree(G)

    fG = graphs.G_to_friendships_graph(G, args.spaligner_long_reads_tsv, args.spaligner_db_tsv)
    graphs.filter_G_by_weight(fG, args.w_name, args.filter)

    # Get feature matrix
    # features_tsv = os.path.join(args.outdir, 'features.tsv')
    # X = graphs.get_X(G.nodes, features_tsv)

    persona_graph, persona_id_mapping = CreatePersonaGraph(fG, local_clustering_fn, args.w_name)

    # graphs drawing
    graphs_outdir = os.path.join(args.outdir, 'graphs_out')
    if not os.path.exists(graphs_outdir):
        os.mkdir(graphs_outdir)
    # graphs.plot_graph_components(G, args.w_name, graphs_outdir, n=4)
    # graphs.plot_graph_components(fG, args.w_name, graphs_outdir, n=4)
    # graphs.plot_graph_components(persona_graph, args.w_name, graphs_outdir, n=10)

    non_overlapping_clustering = list(global_clustering_fn(persona_graph, weight=args.w_name))
    # evaluating_clustering.plot_graph_clusters(persona_graph, non_overlapping_clustering, graphs_outdir)
    # evaluating_clustering.plot_components_clusters(persona_graph, non_overlapping_clustering, args.w_name, graphs_outdir, n=100)

    clustering = PersonaOverlappingClustering(non_overlapping_clustering, persona_id_mapping, 1)

    p_clustering_tsv = os.path.join(args.outdir, 'persona_clustering.tsv')
    evaluating_clustering.write_clustering(non_overlapping_clustering, p_clustering_tsv)

    clustering_tsv = os.path.join(args.outdir, 'clustering.tsv')
    evaluating_clustering.write_clustering(clustering, clustering_tsv)

    nx.write_edgelist(persona_graph, os.path.join(args.outdir, 'persona_graph.tsv'))

    persona_to_node_tsv = os.path.join(args.outdir, 'persona_graph_mapping.tsv')
    with open(persona_to_node_tsv, 'w') as outfile:
        for persona_node, original_node in persona_id_mapping.items():
            outfile.write('{} {}\n'.format(persona_node, original_node))

    '''print('Embedding...')
    embedding = do_embedding(fG, persona_graph, persona_id_mapping,
                             embedding_dim=16, walk_length=10, num_walks_node=40,
                             constraint_learning_rate_scaling_factor=0.1, iterations=10,
                             seed=42)

    # output embedding
    emb_outdir = os.path.join(outdir, 'embedding_out')
    if not os.path.exists(emb_outdir):
        os.mkdir(emb_outdir)

    p_emb_tsv = os.path.join(emb_outdir, 'persona_embedding.tsv')
    embedding['persona_model'].save_word2vec_format(open(p_emb_tsv, 'wb'))
    p_emb_tsv = remove_regular_model(p_emb_tsv, os.path.join(emb_outdir, 'persona_embedding.clear.tsv'))

    # optional output
    embedding['regular_model'].save_word2vec_format(open(os.path.join(emb_outdir, 'embedding_prior.tsv'), 'wb'))

    tot_emb_df = get_total_emb(p_emb_tsv, features_tsv, persona_to_node_tsv)

    visualising_embedding.visualize_embedding(tot_emb_df, persona_to_node_tsv,
                                              spaligner_ground_truth_tsv, p_clustering_tsv,
                                              gfa, fG, emb_outdir)'''

    ground_truth_clustering_tsv = \
        spaligner_parser.spaligner_to_clustering_tsv(args.spaligner_ground_truth_tsv,
                                                     os.path.join(args.outdir, 'ground_truth_clustering.tsv'),
                                                     fG)
    evaluating_clustering.evaluate_clustering(clustering_tsv, ground_truth_clustering_tsv, args.outdir)


def run_with_cProfile():
    import cProfile

    import pstats

    cProfile.run('main()', os.path.join('ml.cprofile'))
    p = pstats.Stats('ml.cprofile')
    p.strip_dirs().sort_stats('cumulative').print_stats(10)


if __name__ == '__main__':
    import random
    random.seed(42)
    import numpy as np
    np.random.seed(42)

    start_time = time.time()

    run_with_cProfile()

    print("Cumtime: %.3f h" % ((time.time() - start_time) / 60 / 60))

