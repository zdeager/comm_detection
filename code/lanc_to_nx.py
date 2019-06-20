#import lancichinetti network/community data into NetworkX 
#
# author: zde

import networkx as nx
# converts network.dat into NetworkX graph
def lanc_to_nx(filename):
	g = nx.Graph()
	network_data = open(filename, 'r')
	for line in network_data:
		edge = line.split()
		g.add_edge(edge[0],edge[1])
	return g

# create true division dictionary using community.dat
def lanc_to_group(filename):
	communities = {}

	comm_data = open(filename, 'r')
	for line in comm_data:
		comm = line.split()
		if (int(comm[1]) not in communities.keys()):
			communities[int(comm[1])] = []
			communities[int(comm[1])].append(comm[0])
		else:
			communities[int(comm[1])].append(comm[0])

	return communities

# compares the true division of community.dat file to that found using
# various comm. detection algorithms
def compare_division(true_divisions, divisions):
	found_divisions  = {}
	for division in divisions:
		max_match = 0, None
		# Find which true division this sub division best matches
		for true_division_key in xrange(1, len(true_divisions.keys()) + 1, 1):
			tmp = len(set(division).intersection(set(true_divisions.get(true_division_key))))
			if (tmp > max_match[0]):
				max_match = tmp, true_division_key

		if (max_match[1] not in found_divisions.keys()):
			found_divisions[max_match[1]] = []
			found_divisions[max_match[1]].extend(division)
		else:
			found_divisions[max_match[1]].extend(division)
	return found_divisions

# computes the number of misplaced nodes based on compare_division()
def incorrect_nodes(true_division, found_division):
	incorrect_nodes = set()

	for key in true_division.keys():
		if (found_division.get(key) == None):
			set_diff = set(true_division.get(key))
		else:
			set_diff = (set(true_division.get(key)).union(set(found_division.get(key))) - \
					set(true_division.get(key)).intersection(set(found_division.get(key))))
		for node in set_diff:
			incorrect_nodes.add(node)
	return incorrect_nodes