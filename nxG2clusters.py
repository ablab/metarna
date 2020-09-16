# nxG2clusters.py assembly_graph_with_scaffolds.gfa outdir

import sys
import os

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

import gfa_parser
import spaligner_parser
import visualising_embedding
import evaluating_clustering


local_clustering_fn = _CLUSTERING_FN['modularity']
global_clustering_fn = _CLUSTERING_FN['modularity']

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

def plot_graph_components(G, outdir, name='G', n=4):
    options = {'with_labels': True,
               'pos': nx.spring_layout(G),
               'font_size': 5}
    largest_components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[:n]
    for i, component in enumerate(largest_components):
        nx.draw(component, **options)
        plt.savefig(os.path.join(outdir, '{}.component_{}.png'.format(name, i)))
        plt.clf()


def main():
    gfa = sys.argv[1]
    spaligner_ground_truth_tsv = sys.argv[2]
    spaligner_long_reads_tsv = sys.argv[3]
    k = int(sys.argv[4])
    outdir = sys.argv[5]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    G = gfa_parser.gfa_to_G(gfa, k)

    # G = get_tst_G(G)

    # Get feature matrix
    features_tsv = os.path.join(outdir, 'features.tsv')
    X = gfa_parser.get_X(G.nodes, features_tsv)

    fG = G.to_undirected()
    plot_graph_components(fG, outdir, n=4)

    # fG.add_edges_from(gfa_parser.get_friendships(G))

    fG.add_edges_from(gfa_parser.get_friendships_from_long_reads(spaligner_long_reads_tsv, fG))
    plot_graph_components(fG, outdir, name='fG', n=4)

    persona_graph, persona_id_mapping = CreatePersonaGraph(fG, local_clustering_fn)
    plot_graph_components(persona_graph, outdir, name='persona', n=10)

    non_overlapping_clustering = list(global_clustering_fn(persona_graph))

    clustering = PersonaOverlappingClustering(non_overlapping_clustering, persona_id_mapping, 2)

    p_clustering_tsv = os.path.join(outdir, 'persona_clustering.tsv')
    with open(p_clustering_tsv, 'w') as outfile:
        for cluster in non_overlapping_clustering:
            outfile.write(' '.join([str(x) for x in cluster]) + '\n')

    clustering_tsv = os.path.join(outdir, 'clustering.tsv')
    with open(clustering_tsv, 'w') as outfile:
        for cluster in clustering:
            outfile.write(' '.join([str(x) for x in cluster]) + '\n')

    nx.write_edgelist(persona_graph, os.path.join(outdir, 'persona_graph.tsv'))

    persona_to_node_tsv = os.path.join(outdir, 'persona_graph_mapping.tsv')
    with open(persona_to_node_tsv, 'w') as outfile:
        for persona_node, original_node in persona_id_mapping.items():
            outfile.write('{} {}\n'.format(persona_node, original_node))

    print('Embedding...')
    embedding = do_embedding(fG, persona_graph, persona_id_mapping,
                             embedding_dim=16, walk_length=10, num_walks_node=40,
                             constraint_learning_rate_scaling_factor=0.1, iterations=10,
                             seed=42)

    # output embedding
    p_emb_tsv = os.path.join(outdir, 'persona_embedding.tsv')
    embedding['persona_model'].save_word2vec_format(open(p_emb_tsv, 'wb'))
    p_emb_tsv = remove_regular_model(p_emb_tsv, os.path.join(outdir, 'persona_embedding.clear.tsv'))

    # optional output
    embedding['regular_model'].save_word2vec_format(open(os.path.join(outdir, 'embedding_prior.tsv'), 'wb'))

    tot_emb_df = get_total_emb(p_emb_tsv, features_tsv, persona_to_node_tsv)

    visualising_embedding.visualize_embedding(tot_emb_df, persona_to_node_tsv,
                                              spaligner_ground_truth_tsv, p_clustering_tsv,
                                              gfa, fG, outdir)

    ground_truth_clustering_tsv = os.path.join(outdir, 'ground_truth_clustering.tsv')
    spaligner_parser.spaligner_to_clustering_tsv(spaligner_ground_truth_tsv, ground_truth_clustering_tsv, fG)
    evaluating_clustering.evaluate_clustering(clustering_tsv, ground_truth_clustering_tsv)


if __name__ == '__main__':
    import random
    random.seed(42)
    import numpy as np
    np.random.seed(42)

    main()
