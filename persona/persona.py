# coding=utf-8
#
# Modifications copyright (C) 2020 Saint Petersburg State University
#
# Copyright 2020 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

r"""Implementation of the Ego-splitting Clustering Framework.

===============================

This is part of the implementation accompanying the WWW 2019 paper, [_Is a
Single Embedding Enough? Learning Node Representations that Capture Multiple
Social Contexts_](https://ai.google/research/pubs/pub46238).

The code in this file allows to create persona graphs, and to obtain overlapping
clusters using the persona graph method defined in the KDD 2017 paper
[_Ego-splitting Framework: from Non-Overlapping to Overlapping
Clusters_](http://epasto.org/papers/kdd2017.pdf).

Citing
------
If you find _Persona Embedding_ useful in your research, we ask that you cite
the following paper:
> Epasto, A., Perozzi, B., (2019).
> Is a Single Embedding Enough? Learning Node Representations that Capture
Multiple Social Contexts.
> In _The Web Conference_.

Example execution
------
python3 -m graph_embedding.persona.persona
  --input_graph=${graph} \
  --output_clustering=${clustering_output}

Where ${graph} is the path to a text file containing the graph and
${clustering_output} is the path to the output clustering.

The graph input format is a text file containing one edge per row represented
as its pair of node ids. The graph is supposed to be directed.
For instance the file:
1 2
2 3

The output clustering format is a text file containing for each row one
(overlapping) cluster represented as the space-separted list of node ids in the
cluster.

The code uses two clustering algorithms local_clustering_method and
global_clustering_method respectively in each egonet and to split the persona ]
graph. The follow three options are allowed at the moment:
connected_components: the standard connected component algorithm.
label_prop: a label propagation based algorithm
            (nx.label_prop.label_propagation_communities).
modularity: an algorithm optimizing modularity
            (nx.modularity.greedy_modularity_communities).
"""

import collections
import itertools
from absl import app
from absl import flags
from boltons.queueutils import HeapPriorityQueue
import networkx as nx
import networkx.algorithms.community.label_propagation as label_prop
import networkx.algorithms.community.modularity_max as modularity
import networkx.algorithms.components.connected as components


_CLUSTERING_FN = {
    'label_prop': label_prop.label_propagation_communities,
    'modularity': modularity.greedy_modularity_communities,
    'connected_components': components.connected_components
}

flags.DEFINE_string(
    'input_graph', None,
    'The input graph path as a text file containing one edge per row, as the '
    'two node ids u v of the edge separated by a whitespace. The graph is '
    'assumed to be directed. For example the file:\n1 2\n2 3\n')

flags.DEFINE_string(
    'output_clustering', None,
    'output path for the overlapping clustering. The clustering is output as a '
    'text file where each row is a cluster, represented as the space-separated '
    'list of its node ids.')

flags.DEFINE_enum(
    'local_clustering_method', 'label_prop', _CLUSTERING_FN.keys(),
    'The method used for clustering the egonets of the graph. The options are '
    '"label_prop", "modularity" or "connected_components".')

flags.DEFINE_enum(
    'global_clustering_method', 'label_prop', _CLUSTERING_FN.keys(),
    'The method used for clustering the persona graph. The options are '
    'label_prop, modularity or connected_components.')

flags.DEFINE_integer('min_cluster_size', 5,
                     'Minimum size for an overlapping cluster to be output.')

flags.DEFINE_string(
    'output_persona_graph', None,
    'If set, it outputs the persona graph in the same format of the input graph'
    ' (text file).')
flags.DEFINE_string(
    'output_persona_graph_mapping', None,
    'If set, outputs the mapping of persona graphs ids to original graph '
    'ids as a text file where each row represents a node and it has two '
    'space-separated columns. The first column is the persona node id, while '
    'the second column is the original node id')

FLAGS = flags.FLAGS


def node_neighbor_to_persona_id(node, egonet, persona_id_counter, clustering_fn,
                                node_neighbor_to_persona_id_map, persona_to_original_mapping):
    partitioning = clustering_fn(egonet)  # Clustering the egonet.
    seen_neighbors = set()
    # Process each of the egonet's local clusters.
    for partition in partitioning:
        persona_id = next(persona_id_counter)
        persona_to_original_mapping[persona_id] = node
        for neighbor in partition:
            node_neighbor_to_persona_id_map[node][neighbor] = persona_id
            assert neighbor not in seen_neighbors
            seen_neighbors.add(neighbor)


def CreatePersonaGraph(graph,
                       clustering_fn=modularity.greedy_modularity_communities,
                       persona_start_id=0):
  """The function creates the persona graph.

  Args:
    graph: Directed graph represented as a dictionary of lists that maps each
      node id its list of neighbor ids;
    clustering_fn: A non-overlapping clustering algorithm function that takes in
      input a nx.Graph and outputs the a clustering. The output format is a list
      containing each partition as element. Each partition is in turn
      represented as a list of node ids. The default function is the networkx
      label_propagation_communities clustering algorithm.
    persona_start_id: The starting id (int) to use for the persona id

  Returns:
    A pair of (graph, mapping) where "graph" is an nx.Graph instance of the
    persona graph (which contains different nodes from the original graph) and
    "mapping" is a dict of the new node ids to the node ids in the original
    graph.The persona graph as nx.Graph, and the mapping of persona nodes to
    original node ids.
  """
  in_egonets = CreateEgonets(graph, direction='in')
  out_egonets = CreateEgonets(graph, direction='out')

  persona_graph = nx.Graph()
  persona_to_original_mapping = dict()
  successor_persona_id_map = collections.defaultdict(dict)
  predecessor_persona_id_map = collections.defaultdict(dict)

  # Next id to allacate in persona graph.
  persona_id_counter = itertools.count(start=persona_start_id)

  for node in graph.nodes():
      # Separate clustering the egonet in forward direction (for output node edges)
      node_neighbor_to_persona_id(node, out_egonets[node], persona_id_counter, clustering_fn,
                                  successor_persona_id_map, persona_to_original_mapping)
      # And separate clustering for input node edges (backward direction)
      node_neighbor_to_persona_id(node, in_egonets[node], persona_id_counter, clustering_fn,
                                  predecessor_persona_id_map, persona_to_original_mapping)

  for u in graph.nodes():  # Process mapping to create persona graph.
    for v in graph.successors(u):
      if v == u:
        continue
      # Since v is successor for u...
      assert v in successor_persona_id_map[u]
      u_p = successor_persona_id_map[u][v]
      # ... u is predecessor for v
      assert u in predecessor_persona_id_map[v]
      v_p = predecessor_persona_id_map[v][u]
      persona_graph.add_edge(u_p, v_p)

  return persona_graph, persona_to_original_mapping


def CreateEgonets(graph, direction):
  """Given a graph, construct all the egonets of the graph.

  Args:
    graph: a nx.Graph instance for which the egonets have to be constructed.

  Returns:
    A dict mapping each node id to an instance of nx.Graph which represents the
    egonet for that node.
  """
  assert direction in ['in', 'out']
  if direction == 'out':
      nbrs = graph.successors
  else:
      nbrs = graph.predecessors

  ego_egonet_map = collections.defaultdict(nx.Graph)

  edge_set = set(graph.edges)

  for node in graph.nodes:
    for neighbor in nbrs(node):
      if neighbor == node:
        continue
      ego_egonet_map[node].add_node(neighbor)

    for u in nbrs(node):
      for v in nbrs(node):
        if (u, v) in edge_set:
          ego_egonet_map[node].add_edge(u, v)

  return ego_egonet_map


def PersonaOverlappingClustering(graph, local_clustering_fn,
                                 global_clustering_fn, min_component_size):
  """Computes an overlapping clustering of graph using the Ego-Splitting method.

  Args:
    graph: a networkx graph for which the egonets have to be constructed.
    local_clustering_fn: method used for clustering the egonets.
    global_clustering_fn: method used for clustering the persona graph.
    min_component_size: minimum size of a cluster to be output.

  Returns:
    The a overlapping clustering (list of sets of node ids), the persona graph
    (nx.Graph) and the persona node
    id mapping (dictionary of int to string) .
  """

  persona_graph, persona_id_mapping = CreatePersonaGraph(graph, local_clustering_fn)
  non_overlapping_clustering = global_clustering_fn(persona_graph)
  overlapping_clustering = set()
  for cluster in non_overlapping_clustering:
    if len(cluster) < min_component_size:
      continue
    cluster_original_graph = set([persona_id_mapping[c] for c in cluster])
    cluster_original_graph = list(cluster_original_graph)
    cluster_original_graph.sort()
    overlapping_clustering.add(tuple(cluster_original_graph))
  return list(overlapping_clustering), persona_graph, persona_id_mapping


def main(argv=()):
  del argv  # Unused.
  graph = nx.read_edgelist(FLAGS.input_graph, create_using=nx.Graph)

  local_clustering_fn = _CLUSTERING_FN[FLAGS.local_clustering_method]
  global_clustering_fn = _CLUSTERING_FN[FLAGS.global_clustering_method]

  clustering, persona_graph, persona_id_mapping = PersonaOverlappingClustering(
      graph, local_clustering_fn, global_clustering_fn, FLAGS.min_cluster_size)

  with open(FLAGS.output_clustering, 'w') as outfile:
    for cluster in clustering:
      outfile.write(' '.join([str(x) for x in cluster]) + '\n')

  if FLAGS.output_persona_graph is not None:
    nx.write_edgelist(persona_graph, FLAGS.output_persona_graph)
  if FLAGS.output_persona_graph_mapping is not None:
    with open(FLAGS.output_persona_graph_mapping, 'w') as outfile:
      for persona_node, original_node in persona_id_mapping.items():
        outfile.write('{} {}\n'.format(persona_node, original_node))
  return 0


if __name__ == '__main__':
  flags.mark_flag_as_required('input_graph')
  flags.mark_flag_as_required('output_clustering')
  app.run(main)
