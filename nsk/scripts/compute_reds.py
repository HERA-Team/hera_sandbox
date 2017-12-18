import numpy as np
import functools

def compute_reds(bls, antpos, tol=2.0):
	"""
	compute redundant baselines

	Parameters:
	-----------
	bls : baseline list, list of antenna pair tuples

	antpos : dictionary, antennas integers as keys, baseline vectors as values

	tol : float, tolerance for redundant baseline determination

	Output:
	-------
	red_bls : redundant baseline list of input bls list
		ordered by smallest separation to longest separation	
	"""
	red_bl_vecs = []
	red_bl_dists = []
	red_bls = []
	for i, k in enumerate(bls):
		bl_vec = antpos[k[1]] - antpos[k[0]]
		unique = map(lambda x: np.linalg.norm(bl_vec - x) > tol, red_bl_vecs)
		if len(unique) == 0 or functools.reduce(lambda x, y: x*y, unique) == 1:
			red_bl_vecs.append(bl_vec)
			red_bl_dists.append(np.linalg.norm(bl_vec))
			red_bls.append([k])
		else:
			red_id = np.where(np.array(unique) == False)[0][0]
			red_bls[red_id].append(k)

	red_bls = list(np.array(red_bls)[np.argsort(red_bl_dists)])
	return red_bls


