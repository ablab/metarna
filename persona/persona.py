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
as its pair of node ids. The graph is supposed to be undirected.
For instance the file:
1 2
2 3
represents the triangle 1, 2, 3.

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
# import networkx.algorithms.community.label_propagation as label_prop

import graphs


def CreatePersonaGraph(graph, clustering_fn, weight_name, persona_start_id=0):
  """The function creates the persona graph.

  Args:
    graph: Undirected graph represented as a dictionary of lists that maps each
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
  egonets = CreateEgonets(graph)
  node_neighbor_persona_id_map = collections.defaultdict(dict)
  persona_graph = nx.OrderedGraph(name='persona')
  persona_to_original_mapping = dict()

  # Next id to allacate in persona graph.
  # persona_id_counter = itertools.count(start=persona_start_id)

  for u, egonet in egonets.items():
    if not egonet.edges:
      partitioning = [[node] for node in egonet.nodes]
    else:
      partitioning = clustering_fn(egonet, weight=weight_name)  # Clustering the egonet.
    seen_neighbors = set()  # Process each of the egonet's local clusters.
    for p_num, partition in enumerate(partitioning):
      # persona_id = next(persona_id_counter)
      persona_id = u + '_' + str(p_num)
      persona_to_original_mapping[persona_id] = u
      for v in partition:
        node_neighbor_persona_id_map[u][v] = persona_id
        assert v not in seen_neighbors
        seen_neighbors.add(v)
  for u in graph.nodes():  # Process mapping to create persona graph.
    for v in graph.neighbors(u):
      if v == u:
        continue
      assert v in node_neighbor_persona_id_map[u]
      u_p = node_neighbor_persona_id_map[u][v]
      assert u in node_neighbor_persona_id_map[v]
      v_p = node_neighbor_persona_id_map[v][u]
      personal_edge_data = get_personal_edge_data(graph, u, v, node_neighbor_persona_id_map, weight_name)
      persona_graph.add_edge(u_p, v_p, **personal_edge_data)
  # graphs.write_G_statistics(persona_graph)
  return persona_graph, persona_to_original_mapping


def get_personal_edge_data(graph, u, v, node_neighbor_persona_id_map, personalized_weight):
  splitted_edge_data = graph.get_edge_data(u, v).copy()
  num_v = len(node_neighbor_persona_id_map[v])
  num_u = len(node_neighbor_persona_id_map[u])
  if num_u == num_v and num_u >= 2:
    splitted_edge_data[personalized_weight] = splitted_edge_data[personalized_weight] * 1.0 / num_u
  return splitted_edge_data


def CreateEgonets(graph):
  """Given a graph, construct all the egonets of the graph.

  Args:
    graph: a nx.Graph instance for which the egonets have to be constructed.

  Returns:
    A dict mapping each node id to an instance of nx.Graph which represents the
    egonet for that node.
  """

  # This is used to not replicate the work for nodes that have been already
  # analyzed..
  completed_nodes = set()
  ego_egonet_map = collections.defaultdict(nx.OrderedGraph)

  # To reducing the running time the nodes are processed in increasing order of
  # degree.
  degrees_pq = HeapPriorityQueue()
  curr_degree = {}
  for node in graph.nodes:
    degrees_pq.add(node, -graph.degree[node])
    curr_degree[node] = graph.degree[node]

  # Ceating a set of the edges for fast membership testing.
  edge_set = set(graph.edges)

  while degrees_pq:
    node = degrees_pq.pop()
    # Update the priority queue decreasing the degree of the neighbor nodes.
    for neighbor in graph.neighbors(node):
      if neighbor == node:
        continue
      ego_egonet_map[node].add_node(
          neighbor)  # even if it is not part of a triangle it is there.
      # We decrease the degree of the nodes still not processed.
      if neighbor not in completed_nodes:
        curr_degree[neighbor] -= 1
        degrees_pq.remove(neighbor)
        degrees_pq.add(neighbor, -curr_degree[neighbor])

    # Construct egonet of u by enumerating all triangles to which u belong
    # because each edge in a triangle is an edge in the egonets of the triangle
    # vertices  and vice versa.
    not_removed = []
    for neighbor in graph.neighbors(node):
      if neighbor not in completed_nodes:
        not_removed.append(neighbor)
    for pos_u, u in enumerate(not_removed):
      for v in not_removed[pos_u + 1:]:
        if (u, v) in edge_set or (v, u) in edge_set:
          ego_egonet_map[node].add_edge(u, v, **graph.get_edge_data(u, v))
          ego_egonet_map[u].add_edge(node, v, **graph.get_edge_data(node, v))
          ego_egonet_map[v].add_edge(u, node, **graph.get_edge_data(u, node))

    completed_nodes.add(node)
  return ego_egonet_map


def PersonaOverlappingClustering(non_overlapping_clustering, persona_id_mapping, min_component_size):
  """Computes an overlapping clustering of graph using the Ego-Splitting method.

  Args:
    non_overlapping_clustering: persona graph clustering
    persona_id_mapping: a dict of the persona node ids to the node ids in the original
    graph
    min_component_size: minimum size of a cluster to be output.

  Returns:
    The overlapping clustering (list of sets of node ids)
  """
  overlapping_clustering = set()
  for cluster in non_overlapping_clustering:
    if len(cluster) < min_component_size:
      continue
    cluster_original_graph = set([persona_id_mapping[persona] for persona in cluster])
    cluster_original_graph = list(cluster_original_graph)
    cluster_original_graph.sort()
    overlapping_clustering.add(tuple(cluster_original_graph))
  return list(overlapping_clustering)