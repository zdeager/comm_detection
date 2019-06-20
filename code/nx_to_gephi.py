# Convert NetworkX graph to a .csv file that can be opened in Gephi.
#
# Author: zde

import numpy # numpy is required (http://www.numpy.org/)
import networkx as nx

# Convert numpy matrix to .csv file that can be read by Gephi
# Params: G - NetworkX graph
# Outputs: out.csv - csv file that can be read by Gephi (File -> Open -> 'out.csv')
def nx_to_gephi(G, filename):
	adjmat_list = nx.adjacency_matrix(G).tolist()

	file = open("%s.csv" % filename, "w")

	temp = ";"
	for node in G.nodes():
		temp += str(node) + ';'
	temp = temp[:-1]
	temp += '\n'
	file.write(temp)


	for node in G.nodes():
		temp = str(node) + ';'
		for entry in adjmat_list[G.nodes().index(node)]:
			temp += str(int(entry)) + ';'
		temp = temp[:-1]
		temp += '\n'
		file.write(temp)

	file.close()