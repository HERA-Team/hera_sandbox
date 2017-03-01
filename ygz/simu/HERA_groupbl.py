import numpy as np
import itertools
import w_opp
from joblib import Parallel, delayed
"""
This file groups the HERA baselines into equivalency classes
"""
FILE = "../calfiles/HERA_antconfig/antenna_positions_37.dat"
def get_plotsense_dict(file=FILE, NTOP=None):
	antpos = np.loadtxt(file)
	NANTS = antpos.shape[0]

	bl_gps = {}
	for i in xrange(NANTS-1):
		a1pos = antpos[i]
		for j in xrange(i+1,NANTS):
			a2pos = antpos[j]
			blx, bly, blz = a2pos-a1pos
			has_entry = False
			for key in bl_gps.keys():
				if np.hypot(key[0]-blx, key[1]-bly) < 1:
					has_entry = True
					bl_gps[key] = bl_gps.get(key, []) + [(i,j)]
					break
			if not has_entry:
				bl_gps[(blx, bly)] = [(i,j)]

	n_unique = len(bl_gps.keys())
	n_total = NANTS*(NANTS-1)/2
	
	print "Found %d classes among %d total baselines" % (n_unique, n_total)
	print "Sorting dictionary"
	sorted_keys = sorted(bl_gps.keys(), key=lambda x: np.hypot(*x))
	if NTOP is None: 
		NTOP = n_unique
	key_pool = np.arange(NTOP)
	top_dict = {}
	for i, key in enumerate(sorted_keys[:NTOP]):
		mult = len(bl_gps[key])
		label = str(bl_gps[key][0][0])+'_'+str(bl_gps[key][0][1])
		top_dict[label] = ((key[0], key[1], 0.), mult) #dict[example_bl] = ((blx, bly, blz), multiplicity)
	return top_dict, bl_gps

def get_bl_comb(top_dict):
	combs = []
	for sing in itertools.combinations(top_dict.iteritems(), 1):
		combs.append((sing[0], sing[0]))
	for pair in itertools.combinations(top_dict.iteritems(), 2):
		combs.append(pair)
	return combs



def run_opp(i, comb, outfile):
	#comb=(('40', (44, 26)), ('41', (44, 38)))
	file = open(outfile, 'a')
	label1, label2 = comb[0][0], comb[1][0]
	bl1coords, bl2coords = comb[0][1][0], comb[1][1][0]
	multiplicity = comb[0][1][1]*comb[1][1][1]

	peak,dT = WS.w_opp(bl1coords=bl1coords,bl2coords=bl2coords)
	line = ', '.join([str(i), label1,label2,str(dT),str(np.abs(peak)),str(multiplicity)])
	print line
	file.write(line+'\n')
	file.close()



if __name__=="__main__":
	top_dict, blgps = get_plotsense_dict()
	combs = get_bl_comb(top_dict)
	CAL = 'psa6622_v003'
	WS = w_opp.OppSolver(fq=.15, cal=CAL)
	OUT = 'HERA_37_opp.csv'
	NJOBS = 4

	print 'Starting Opp with %d instances on %d jobs' % (len(combs), NJOBS)
	print ',sep,sep2,dT,peak,mult'
	Parallel(n_jobs=4)(delayed(run_opp)(i, comb, OUT) for i, comb in enumerate(combs))

	import IPython; IPython.embed()