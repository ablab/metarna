# nxG2clusters.py assembly_graph_with_scaffolds.gfa outdir

import sys
import os

import networkx as nx

from persona.persona import CreatePersonaGraph, PersonaOverlappingClustering
from persona.persona import _CLUSTERING_FN
from persona.splitter import do_embedding

import gfa2nxG

local_clustering_fn = _CLUSTERING_FN['weakly_connected_components']
global_clustering_fn = _CLUSTERING_FN['weakly_connected_components']

gfa = sys.argv[1]
G = gfa2nxG.gfa_to_G(gfa)

outdir = sys.argv[2]

# G for testing
# G = nx.DiGraph()
# G.add_nodes_from(['a', 'b', 'c', 'd', 'e'])
# G.add_edges_from([('a', 'b'), ('a', 'c'), ('b', 'd'), ('b', 'e'), ('c', 'd'), ('c', 'e')])

# path1 930004-,278546-,36185+,278990+,283130+,352975-,37703+
# path2 930004-,239212-,36185+,365256-,283130+,352975-,37703+
# nodes_tst = ['36185+', '37703+', '239212-', '278546-', '278990+',
#              '283130+', '352975-', '365256-', '930004-', '2326645-']
# G_tst = G.subgraph(nodes_tst).copy()
# G_tst.add_edges_from([('36185+', '278990+'), ('283130+', '352975-'), ('36185+', '365256-'), ('283130+', '2326645-')])
# X = gfa2nxG.get_X(G_tst)
# A = gfa2nxG.get_A(G_tst)

persona_graph, persona_id_mapping = CreatePersonaGraph(G, local_clustering_fn)
non_overlapping_clustering = list(global_clustering_fn(persona_graph))
clustering = PersonaOverlappingClustering(non_overlapping_clustering, persona_id_mapping, 0)

with open(os.path.join(outdir, 'persona_clustering.tsv'), 'w') as outfile:
    for cluster in non_overlapping_clustering:
        outfile.write(' '.join([str(x) for x in cluster]) + '\n')

with open(os.path.join(outdir, 'clustering.tsv'), 'w') as outfile:
    for cluster in clustering:
        outfile.write(' '.join([str(x) for x in cluster]) + '\n')

nx.write_edgelist(persona_graph, os.path.join(outdir, 'persona_graph.tsv'))

with open(os.path.join(outdir, 'persona_graph_mapping.tsv'), 'w') as outfile:
    for persona_node, original_node in persona_id_mapping.items():
        outfile.write('{} {}\n'.format(persona_node, original_node))

print('Running splitter...')
embedding = do_embedding(G, persona_graph, persona_id_mapping,
                         embedding_dim=64, walk_length=10, num_walks_node=40,
                         constraint_learning_rate_scaling_factor=0.1, iterations=10,
                         seed=1)

# output embeddings
def remove_regular_model(in_path, out_path):
    fout = open(out_path, 'w')

    with open(in_path, 'r') as fin:
        for line in fin:
            node = line.split()[0]
            if not ('+' in node or '-' in node):
                fout.write(line)

    fout.close()

persona_model = os.path.join(outdir, 'persona_embedding.tsv')
embedding['persona_model'].save_word2vec_format(open(persona_model, 'wb'))
remove_regular_model(persona_model, os.path.join(outdir, 'persona_embedding.clear.tsv'))

# optional output
embedding['regular_model'].save_word2vec_format(open(os.path.join(outdir, 'embedding_prior.tsv'), 'wb'))

