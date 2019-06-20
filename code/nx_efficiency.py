# Implements division of NetworkX graph based into communities based
# on eff. modularity algorithm
#
# Author: zde
# Author: Kristen Bales
import networkx as nx
import numpy
import math

# Returns the leading eigenvector of a modularity matrix using built in numpy solver
# Params - mod_mat : numpy array
# Returns - eig_vec : leading eigenvector
def leading_eig_vec(mod_mat):
	eig_vec_and_vals = numpy.linalg.eigh(mod_mat)
	index = numpy.argmax(numpy.real(eig_vec_and_vals[0])) # Get index of leading eigenvalue
	return numpy.real(eig_vec_and_vals[1][:,index])

# Calculate the eff. modularity of a NetworkX graph from list of NetworkX subgraphs
# Params - g : the original NetworkX graph
#		   subgraphs : the list of subgraphs
# Returns - eff. modularity of g (float)
def nx_eff_modularity(g, subgraphs, p):
	if (g.size() == 0):
		return 0
	E = nx_efficiency_matrix(g, p)
	m = nx_m_tilde(g, p)
	n = g.order()
	k_vec = numpy.zeros(n)

	for i in xrange(n):
		e_sum = 0
		for j in xrange(n):
			e_sum += E[i][j]
		k_vec[i] = e_sum #sum of each row

	indices = []
	f_sum = 0
	for subgraph in xrange(0, len(subgraphs), 1):
		for node in subgraphs[subgraph]:
			indices.append(g.nodes().index(node))
		# Loop over nodes i,j in subgraph pair wise
		for i in indices:
			for j in indices:
				f_sum += E[i][j] - ((k_vec[i] * k_vec[j]) / (2.0 * m))
		indices = []
	return (f_sum / (2.0 * m))

# Returns the value of \tilde{m} from NetworkX graph g
# Params - g : NetworkX graph
# Returns - m_tilde : float
def nx_m_tilde(g, p):
	E = nx_efficiency_matrix(g, p)
	n = g.order()
	m_tilde = 0
	for i in xrange(n):
		for j in xrange(n):
			m_tilde += E[i][j]
	m_tilde = m_tilde/2
	return  m_tilde

# Returns the efficiency matrix (modified to 1/(d^p) to penalize longer paths more)
# Params - g : NetworkX graph
#		   p : int (power to raise distance matrix entries to)
# Returns - e_mat : numpy array
def nx_efficiency_matrix(g, p):
	n = g.order()
	d = nx.floyd_warshall_numpy(g).tolist() # fastest NetworkX all-pairs shortest path
	
	e_mat = numpy.zeros((n,n))
	for i in xrange(0, n, 1):
		for j in xrange(i, n, 1):
			if d[i][j] == 0 or d[i][j] == float('Inf'):
				E_ij = 0
			else:
				E_ij = 1 / math.pow(d[i][j],p) # 1/(d^p)
				#E_ij = 1 / math.exp(d[i][j] - 1) # 1/(d^p)
			e_mat[i][j] = E_ij
			e_mat[j][i] = E_ij
	return e_mat

# Compute/Return eff. modularity matrix for a NetworkX graph
# If calculating B_tilde_g for subgraph of g, g_orig is the the original graph and 
# 			B_orig is the eff. modularity matrix of the original graph that g is a subgraph of
# Params - g : NetworkX graph
#		   g_orig : the original graph that g is a subgraph of
# 		   B_orig : the eff. modularity matrix of the original graph that g is a subgraph of
# Returns - B/B_g : the eff. modularity matrix as a numpy array
def B_tilde(g, p, g_orig, B_orig):
	if (g_orig == None and B_orig == None): # Calculating eff. modularity matrix for orignal graph
		E = nx_efficiency_matrix(g, p)
		n = g.order()
		k_vec = numpy.zeros(n)

		m_tilde = 0
		for i in xrange(n):
			e_sum = 0
			for j in xrange(n):
				e_sum += E[i][j]
			k_vec[i] = e_sum #sum of each row
			m_tilde += e_sum #total sum of E matrix
		m_tilde = m_tilde/2

		E_tilde = numpy.zeros((g.order(),g.order()))
		for i in xrange(n):
			for j in xrange(i, n, 1):
				e_ij =  k_vec[i]*k_vec[j]/(2*m_tilde)
				E_tilde[i][j] = e_ij
				E_tilde[j][i] = e_ij
		return E - E_tilde
	else: # Calculating modularity matrix for subgraph of g_orig
		B_g = numpy.zeros((g.order(),g.order()))
		indices = []
		for node in g:
			indices.append(g_orig.nodes().index(node))
		for i in indices:
			for j in indices:
				if (i != j):
					B_g[indices.index(i)][indices.index(j)] = B_orig[i][j]
				else:
					sum_B_ik = 0
					for k in indices:
						sum_B_ik += B_orig[i][k]
					B_g[indices.index(i)][indices.index(j)] = B_orig[i][j] - sum_B_ik
		return B_g

# Partition a NetworkX graph (if divisible) into communities based on the eff. modularity
# Params - g : NetworkX graph
# 		   n : number of leading eigenvectors to test for optimal partitions
#		   depth : the desired number or partitions to be preformed (limit of recursion)
#		   mod_mat_orig : the modularity matrix of the graph being partitioned
#		   g_orig : the original graph being partitioned
#		   curr_depth : current recursion depth
# Returns - communities : list of subgraphs of g
def nx_divide_graph(g, depth = float("inf"), swap = True, p = 1, mod_mat_orig = None, g_orig = None, curr_depth = 0):
	subgraphs = []

	# Graph has no edges (g is indivisble) or desired number of divisions was reached
	if (g.size() == 0 or depth <= curr_depth): 
		subgraphs.append(g)
		return subgraphs
	# If this is the first iteration, set mod_mat_orig and g_orig
	if (mod_mat_orig == None and g_orig == None):
		mod_mat_orig = B_tilde(g, p, g_orig, mod_mat_orig)
		g_orig = g
		mod_mat = mod_mat_orig
	else: # Else calculate B_tilde_g
		mod_mat = B_tilde(g, p, g_orig, mod_mat_orig)

	# Calculate leading eigenvector
	eig_vec = leading_eig_vec(mod_mat[:])
	
	s = [] # List of s_i values (+/-1)
	for x in xrange(0, len(eig_vec), 1):
		if eig_vec[x] >= 0: # If entry of eigenvector is positive, place node in first community
			s.append(1)
		else: # Else place node in second community
			s.append(-1)

	m = nx_m_tilde(g, p)
	delta_q = (1.0 / (4 * m)) * numpy.dot(numpy.dot(numpy.transpose(s), mod_mat),s) # Calc delta Q
	print(delta_q)
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
	if (best_delta_q > -0.0001 and len(communities[0]) > 0 and len(communities[1]) > 0):
		subgraphs.extend(nx_divide_graph(g.subgraph(communities[0]), depth, swap, p, \
													mod_mat_orig, g_orig, curr_depth + 1))
		subgraphs.extend(nx_divide_graph(g.subgraph(communities[1]), depth, swap, p, \
													mod_mat_orig, g_orig, curr_depth + 1))
	else: # Else, g is indivisble
		subgraphs.append(g)

	return subgraphs