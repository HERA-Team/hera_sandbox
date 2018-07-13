import numpy as np
import functools

def compute_reds(antpos, ex_ants=[], tol=1.0, pol=None):
    """
    compute redundant baselines groups.

    Parameters:
    -----------
    antpos : dictionary, antennas integers as keys, baseline vectors as values

    ex_ants : list of flagged (excluded) antennas

    tol : float, tolerance for redundant baseline determination in units of the baseline vector units

    pols : type=str, polarization string to append to keys

    Output:
    -------
    red_bls : redundant baseline list of input bls list
        ordered by smallest separation to longest separation    
    """
    if type(antpos) is not dict and type(antpos) is not odict:
        raise AttributeError("antpos is not a dictionary type")

    # calculate all permutations
    bls = map(lambda bl: tuple(sorted(bl)), sorted(itertools.combinations(antpos.keys(), 2)))

    # add pol string
    if pol is not None:
        bls = map(lambda bl: bl + tuple([pol]), bls)

    red_bl_vecs = []
    red_bl_dists = []
    red_bls = []
    for i, k in enumerate(bls):
        if k[0] in ex_ants or k[1] in ex_ants:
            continue
        try:
            bl_vec = antpos[k[1]] - antpos[k[0]]
        except KeyError:
            continue
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


