# nxG2clusters.py assembly_graph_with_scaffolds.gfa outdir

import sys
import os

import networkx as nx

from persona.persona import CreatePersonaGraph, PersonaOverlappingClustering
from persona.persona import _CLUSTERING_FN
from persona.splitter import do_embedding

import gfa2nxG
import visualising_PCA_tSNE

local_clustering_fn = _CLUSTERING_FN['weakly_connected_components']
global_clustering_fn = _CLUSTERING_FN['weakly_connected_components']

def remove_regular_model(in_path, out_path):
    fout = open(out_path, 'w')

    with open(in_path, 'r') as fin:
        for line in fin:
            node = line.split()[0]
            if not ('+' in node or '-' in node):
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

def main():
    gfa = sys.argv[1]
    spaligner_tsv = sys.argv[2]
    outdir = sys.argv[3]

    G = gfa2nxG.gfa_to_G(gfa)

    # G_tst = get_tst_G(G)

    # Get feature matrix
    features_tsv = os.path.join(outdir, 'features.tsv')
    X = gfa2nxG.get_X(G, features_tsv)
    # Set labels for nodes
    node_to_db_tsv = os.path.join(outdir, 'node_to_db.tsv')
    G = gfa2nxG.set_node_labels(G, spaligner_tsv, node_to_db_tsv)

    persona_graph, persona_id_mapping = CreatePersonaGraph(G, local_clustering_fn)

    non_overlapping_clustering = list(global_clustering_fn(persona_graph))

    clustering = PersonaOverlappingClustering(non_overlapping_clustering, persona_id_mapping, 0)

    p_clustering_tsv = os.path.join(outdir, 'persona_clustering.tsv')
    with open(p_clustering_tsv, 'w') as outfile:
        for cluster in non_overlapping_clustering:
            outfile.write(' '.join([str(x) for x in cluster]) + '\n')

    with open(os.path.join(outdir, 'clustering.tsv'), 'w') as outfile:
        for cluster in clustering:
            outfile.write(' '.join([str(x) for x in cluster]) + '\n')

    nx.write_edgelist(persona_graph, os.path.join(outdir, 'persona_graph.tsv'))

    persona_to_node_tsv = os.path.join(outdir, 'persona_graph_mapping.tsv')
    with open(persona_to_node_tsv, 'w') as outfile:
        for persona_node, original_node in persona_id_mapping.items():
            outfile.write('{} {}\n'.format(persona_node, original_node))

    print('Embedding...')
    embedding = do_embedding(G, persona_graph, persona_id_mapping,
                             embedding_dim=64, walk_length=10, num_walks_node=40,
                             constraint_learning_rate_scaling_factor=0.1, iterations=10,
                             seed=1)

    # output embedding
    p_emb_tsv = os.path.join(outdir, 'persona_embedding.tsv')
    embedding['persona_model'].save_word2vec_format(open(p_emb_tsv, 'wb'))
    p_emb_tsv = remove_regular_model(p_emb_tsv, os.path.join(outdir, 'persona_embedding.clear.tsv'))

    # optional output
    embedding['regular_model'].save_word2vec_format(open(os.path.join(outdir, 'embedding_prior.tsv'), 'wb'))

    visualising_PCA_tSNE.visualize_embedding(p_emb_tsv, persona_to_node_tsv, node_to_db_tsv, p_clustering_tsv, outdir)


if __name__ == '__main__':
    main()
