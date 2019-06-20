# Implementation of graph efficiency functions for NetworkX Graphs.
# Implements local and global efficiency functions.
# 
# Author: zde

from __future__ import division
import networkx as nx

# Return the local efficiency of each node in the graph G
# Params: G - NetworkX graph
# Returns: local_efficiency - dict
#       	the keys of the dict are the nodes in the graph G and the 
#			values are local efficiencies of each node
def local_efficiencies(G, weight = None):
    if G.is_directed():
        new_graph = nx.DiGraph
    else:
        new_graph = nx.Graph

    efficiencies = {}
    for node in G:
    	# Create local subgraph
        temp_G = nx.Graph() 
        temp_G.add_nodes_from(G.neighbors(node))
        for neighbor in G.neighbors(node):
            for (n1, n2) in G.edges(neighbor):
                if (n1 in temp_G) and (n2 in temp_G):
                    temp_G.add_edge(n1, n2)
        # add wights to subgraph
        if weight is not None:
            for (n1, n2) in temp_G.edges():
                print G[n1][n2]
                for dic in G[n1][n2].values():
                    if dic.keys()[0] == weight:
                        temp_G[n1][n2][weight] = dic[weight]
                        break
        # Find global efficiency of subgraph
        efficiencies[node] = global_efficiency(temp_G)

    return efficiencies

# Return the average local efficiency of all of the nodes in the graph G (efficiency of subgraphs)
# Params: G - NetworkX graph
# Returns: local_efficiency - float
def local_efficiency(G, weight = None):
    eff = local_efficiencies(G, weight)
    total = sum(eff.values())
    N = len(eff)
    return total/N

# Return the global efficiency of the graph G
# Params: G - NetworkX graph
# Returns: global_efficiency - float
def global_efficiency(G, weight = None):
    N = len(G) # get order of graph

    inv_lengths = []
    for node in G:
        # Get path lengths from this node to every reachable node in G
        if weight is None:
            lengths = nx.single_source_shortest_path_length(G, node)
        else:
            lengths = nx.single_source_dijkstra_path_length(G, node, weight=weight)
        # Find inverse of the path length to each reachable node (excluding this node)
        inv = [1/x for x in lengths.values() if x is not 0]
        inv_lengths.extend(inv)

    return sum(inv_lengths)/(N*(N-1))


# Return a graph where the weight of each edge is 1 over the weight of
# the corresponding edge in graph G
# NOTE: the key in the dictionary is 'inv_weight'
# Params: G - NetworkX graph
# Returns: Ginv - NetworkX graph
def invert_weights(G, weight='weight'):
    Ginv = nx.create_empty_copy(G)
    for (n1, n2) in G.edges():
        for dic in G[n1][n2].values():
            if dic.keys()[0] == weight:
                dist = 1/dic[weight]
                break
        Ginv.add_edge(n1, n2, inv_weight = dist)

    return Ginv