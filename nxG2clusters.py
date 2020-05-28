# nxG2clusters.py assembly_graph_with_scaffolds.gfa outdir

import sys
import os

import networkx as nx

from persona.persona import PersonaOverlappingClustering
from persona.persona import _CLUSTERING_FN

import gfa2nxG


local_clustering_fn = _CLUSTERING_FN['label_prop']
global_clustering_fn = _CLUSTERING_FN['label_prop']

gfa = sys.argv[1]
G = gfa2nxG.gfa_to_G(gfa)

# G for testing
# G = nx.DiGraph()
# G.add_nodes_from(['a', 'b', 'c', 'd', 'e'])
# G.add_edges_from([('a', 'b'), ('a', 'c'), ('b', 'd'), ('b', 'e'), ('c', 'd'), ('c', 'e')])

# path1 930004-,278546-,36185+,278990+,283130+,352975-,37703+
# path2 930004-,239212-,36185+,365256-,283130+,352975-,37703+
nodes_tst = ['36185+', '37703+', '239212-', '278546-', '278990+',
             '283130+', '352975-', '365256-', '930004-', '2326645-']
G_tst = G.subgraph(nodes_tst)
X = gfa2nxG.get_X(G_tst)

clustering, persona_graph, persona_id_mapping = \
    PersonaOverlappingClustering(G_tst, local_clustering_fn, global_clustering_fn, 0)

outdir = sys.argv[2]
with open(os.path.join(outdir, 'clustering.tsv'), 'w') as outfile:
    for cluster in clustering:
        outfile.write(' '.join([str(x) for x in cluster]) + '\n')

nx.write_edgelist(persona_graph, os.path.join(outdir, 'persona_graph.tsv'))

with open(os.path.join(outdir, 'persona_graph_mapping.tsv'), 'w') as outfile:
    for persona_node, original_node in persona_id_mapping.items():
        outfile.write('{} {}\n'.format(persona_node, original_node))