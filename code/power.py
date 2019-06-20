# attempt at fast power method
#
# author: zde
def leading_eig_vec(mod_mat, maxIters = 1000, tolerance = 1e-4):
	v = numpy.random.random_integers(0,1,size = len(mod_mat)) # random non-zero vector for initial guess
	iteration = 0; l_old = 0
	dominant_eig = None

	while (iteration <= maxIters):
		z = numpy.dot(mod_mat,v) # transform mod_mat*v
		v = z/numpy.linalg.norm(z) # normalize
		l = numpy.dot(numpy.dot(mod_mat,v).T,v) / (numpy.linalg.norm(v) ** 2) # calc new lambda
		iteration = iteration + 1
		if (abs((l - l_old)/l) < tolerance):
			dominant_eig = (l, v)
			break
		l_old = l

	if (dominant_eig == None): # No convergence
		return None
	elif (dominant_eig[0] < 0): # Subtract lambda*I to mod_mat and redo power method
		mod_mat -= (numpy.identity(len(mod_mat))*dominant_eig[0])
		return leading_eig_vec(mod_mat)
	else:
		return dominant_eig[1]
'''