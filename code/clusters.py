# Creates an returns a NetworkX graph with complete graphs connected by
# the desired number of egdes. Used for testing community detection.
#

import random as r
import networkx as nx 

# Returns m number of K_n subgraphs connected to every other K_n by e edges
# Params - n : order of complete graphs K_n (int)
#          m : number of K_n communities (int)
#          e : number of edges connecting K_n communities (int)
# Returns - g : NetworkX graph of K_n communities
def clusters(n, m, e=1):
	g = nx.Graph()
	for i in xrange(m):
		h = nx.complete_graph(n)
		h = nx.convert_node_labels_to_integers(h,first_label=(i*n))
		g = nx.union(g,h)
	for i in xrange(m):
		for j in xrange(i+1, m):
			for k in xrange(e):
				A = r.randint(i*n, (i+1)*n-1)
				B = r.randint(j*n, (j+1)*n-1)
				g.add_edge(A,B)
	return g