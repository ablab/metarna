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
as its pair of node ids.

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
from boltons.queueutils import HeapPriorityQueue
import networkx as nx


# Create map of persona id for node neighbors in directed graphs
def neighbor_to_persona_id(node, egonet, persona_id_counter, clustering_fn,
                           neighbor_to_persona_id_map, persona_to_original_mapping):
    partitioning = clustering_fn(egonet)  # Clustering the egonet.
    seen_neighbors = set()
    # Process each of the egonet's local clusters.
    for partition in partitioning:
        persona_id = next(persona_id_counter)
        persona_to_original_mapping[persona_id] = node
        for neighbor in partition:
            neighbor_to_persona_id_map[neighbor] = persona_id
            assert neighbor not in seen_neighbors
            seen_neighbors.add(neighbor)

def CreateDirectedPersonaGraph(graph, clustering_fn, persona_start_id=0):
  """The function creates the directed persona graph.

  Args:
    graph: Directed graph;
    clustering_fn: A non-overlapping clustering algorithm function;
    persona_start_id: The starting int id to use for the persona id;

  Returns:
    The persona graph as nx.DiGraph, and the mapping of persona nodes to
    original node ids.
  """
  in_egonets = CreateDirectedEgonets(graph, direction='in')
  out_egonets = CreateDirectedEgonets(graph, direction='out')

  persona_graph = nx.DiGraph()
  persona_to_original_mapping = dict()
  successor_persona_id_map = collections.defaultdict(dict)
  predecessor_persona_id_map = collections.defaultdict(dict)

  # Next id to allocate in persona graph.
  persona_id_counter = itertools.count(start=persona_start_id)
  for node in graph.nodes():
      # Separate clustering the egonet in forward direction (for output node edges)
      neighbor_to_persona_id(node, out_egonets[node], persona_id_counter, clustering_fn,
                             successor_persona_id_map[node], persona_to_original_mapping)
      # And separate clustering for input node edges (backward direction)
      neighbor_to_persona_id(node, in_egonets[node], persona_id_counter, clustering_fn,
                             predecessor_persona_id_map[node], persona_to_original_mapping)

  persona_graph.add_nodes_from(persona_to_original_mapping.keys())
  for u in graph.nodes():  # Process mapping to create persona graph.
    for v in graph.successors(u):
      if v == u:
        continue
      # Since v is successor for u...
      assert v in successor_persona_id_map[u]
      successor = successor_persona_id_map[u][v]
      # ... u is predecessor for v
      assert u in predecessor_persona_id_map[v]
      predecessor = predecessor_persona_id_map[v][u]
      persona_graph.add_edge(successor, predecessor)

  return persona_graph, persona_to_original_mapping

def CreateDirectedEgonets(graph, direction):
  """Given a directed graph, construct all the egonets of the graph.

  Args:
    graph: a nx.diGraph instance for which the egonets have to be constructed.

  Returns:
    A dict mapping each node id to an instance of nx.diGraph which represents the
    egonet for that node.
  """
  assert direction in ['in', 'out']
  if direction == 'out':
      nbrs = graph.successors
  else:
      nbrs = graph.predecessors

  ego_egonet_map = collections.defaultdict(nx.DiGraph)

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