# Implements division of NetworkX graph based into communities based
# on edge modularity algorithm presented in http://www.pnas.org/content/103/23/8577.full.
#
# Author: zde
# Author: Kristen Bales
import networkx as nx
import numpy

# Return the n leading eigenvectors of a modularity matrix
# Params - mod_mat : numpy array
# Returns - eig_vec : leading eigenvector
def leading_eig_vec(mod_mat):
	eig_vec_and_vals = numpy.linalg.eig(mod_mat)
	index = numpy.argmax(numpy.real(eig_vec_and_vals[0])) # Get index of leading eigenvalue
	return numpy.real(eig_vec_and_vals[1][:,index])

# Calculate the edge modularity of a NetworkX graph from list of NetworkX subgraphs
# Params - g : the original NetworkX graph
#		   subgraphs : the list of subgraphs
# Returns - modularity of g (float)
def nx_edge_modularity(g, subgraphs):
	if (g.size() == 0):
		return 0
	g_adj_mat = nx.to_numpy_matrix(g).tolist()
	indices = []
	f_sum = 0
	for subgraph in xrange(0, len(subgraphs), 1):
		for node in subgraphs[subgraph]:
			indices.append(g.nodes().index(node))
		# Loop over nodes i,j in subgraph pair wise
		for i in indices:
			for j in indices:
				f_sum += g_adj_mat[i][j] - \
					((g.degree(g.nodes()[i]) * g.degree(g.nodes()[j])) / (2.0 * g.size()))
		indices = []
	return (f_sum / (2.0 * g.size()))

# Compute/Return modularity matrix for a NetworkX graph
# If calculating B_g for subgraph of g, g_orig is the the original graph and 
# 			B_orig is the modularity matrix of the original graph that g is a subgraph of
# Params - g : NetworkX graph
#		   g_orig : the original graph that g is a subgraph of
# 		   B_orig : the modularity matrix of the original graph that g is a subgraph of
# Returns - B/B_g : the modularity matrix as a numpy array
def nx_modularity_matrix(g, g_orig = None, B_orig = None):
	if (g_orig == None and B_orig == None): # Calculating modularity matrix for orignal graph
		n = g.order()
		m = g.size()
		B = nx.to_numpy_matrix(g).tolist()
		deg_list = [] # List containing degrees of nodes
		for node in g.nodes():
			deg_list.append(g.degree(node))
		# Loop over nodes i,j pair wise
		for i in xrange(0, n, 1):
			for j in xrange(i, n, 1):
				# Only calculate B_ij once since B_ji = B_ij
				B_ij = B[i][j] - (deg_list[i] * deg_list[j] / (2.0 * m))
				B[i][j] = B_ij
				B[j][i] = B_ij
		return B
	else: # Calculating modularity matrix for subgraph of g_orig
		B_g = numpy.zeros((g.order(),g.order()))
		# Get indices of nodes in original graph
		indices = []
		for node in g:
			indices.append(g_orig.nodes().index(node))
		# Loop over nodes i,j in subgraph pair wise
		for i in indices:
			for j in indices:
				if (i != j): # If i != j Kronecker delta is 0
					B_g[indices.index(i)][indices.index(j)] = B_orig[i][j]
				else: # Else Kronecker delta is 1, so calculate sum(B_ik)
					sum_B_ik = 0
					for k in indices:
						sum_B_ik += B_orig[i][k]
					B_g[indices.index(i)][indices.index(j)] = B_orig[i][j] - sum_B_ik
		return B_g

# Partition a NetworkX graph (if divisible) into communities based on the algorithm
# presented in http://www.pnas.org/content/103/23/8577.full
# Params - g : NetworkX graph
# 		   n : number of leading eigenvectors to test for optimal partitions
#		   depth : the desired number or partitions to be preformed (limit of recursion)
#		   mod_mat_orig : the modularity matrix of the graph being partitioned
#		   g_orig : the original graph being partitioned
#		   curr_depth : current recursion depth
# Returns - communities : list of subgraphs of g
def nx_divide_graph(g, depth = float("inf"), swap = True, mod_mat_orig = None, g_orig = None, curr_depth = 0):
	subgraphs = []
	m = g.size()

	# Graph has no edges (g is indivisble) or desired number of divisions was reached
	if (m == 0 or depth <= curr_depth): 
		subgraphs.append(g)
		return subgraphs
	# If this is the first iteration, set mod_mat_orig and g_orig
	if (mod_mat_orig == None and g_orig == None):
		mod_mat_orig = nx_modularity_matrix(g)
		g_orig = g
		mod_mat = mod_mat_orig
	else: # Else calculate B_g
		mod_mat = nx_modularity_matrix(g, g_orig, mod_mat_orig)

	# Calculate leading eigenvector
	eig_vec = leading_eig_vec(mod_mat[:])
	s = [] # List of s_i values (+/-1)
	for x in xrange(0, len(eig_vec), 1):
		if eig_vec[x] >= 0: # If entry of eigenvector is positive, place node in first community
			s.append(1)
		else: # Else place node in second community
			s.append(-1)

	delta_q = (1.0 / (4 * m)) * numpy.dot(numpy.dot(numpy.transpose(s), mod_mat),s) # Calc delta Q
	
	# Perform node swapping
	best_delta_q = delta_q
	if (swap):
		swapped_nodes = []
		while (True):
			temp_configs = {}
			for index in xrange(0, len(s), 1): # Swap each node and calcluate delta q
				if (index not in swapped_nodes): # Only swap node once
					s[index] *= -1 # flip s[index]
					delta_q = (1.0 / (4 * m)) * numpy.dot(numpy.dot(numpy.transpose(s), mod_mat),s) # Calc delta Q
					temp_configs[delta_q] = index
					s[index] *= -1 # flip back
			temp_max_q = max(temp_configs.keys())
			if (temp_max_q > best_delta_q): # Keep node swap which increased delta q the most
				s[temp_configs[temp_max_q]] *= -1
				swapped_nodes.append(temp_configs[temp_max_q])
				best_delta_q = temp_max_q
			else: # No increase in delta q is possible
				break;
	
	# Finalize the node division
	communities = [[],[]]
	for x in xrange(0, len(s), 1):
		if s[x] >= 0: # If entry of eigenvector is positive, place node in first community
			communities[0].append(g.nodes()[x])
		else: # Else place node in second community
			communities[1].append(g.nodes()[x])

	# Recursively divide the optimal partition if delta Q is positive
	if (best_delta_q > 0 and len(communities[0]) > 0 and len(communities[1]) > 0):
		subgraphs.extend(nx_divide_graph(g.subgraph(communities[0]), depth, swap, \
													mod_mat_orig, g_orig, curr_depth + 1))
		subgraphs.extend(nx_divide_graph(g.subgraph(communities[1]), depth, swap, \
													mod_mat_orig, g_orig, curr_depth + 1))
	else: # Else, g is indivisble
		subgraphs.append(g)

	return subgraphs