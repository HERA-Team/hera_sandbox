import numpy as np, pandas as pd
import itertools, pickle
import w_opp
from joblib import Parallel, delayed
"""
This file groups the HERA baselines into equivalency classes
"""
def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def get_plotsense_dict(file, NTOP=None, NANTS=None):
	antpos = np.loadtxt(file)
	if NANTS is None:
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

def get_bl_comb(top_dict, alpha=None):
	combs = []
	for sing in itertools.combinations(top_dict.iteritems(), 1):
		combs.append((sing[0], sing[0]))
	for pair in itertools.combinations(top_dict.iteritems(), 2):
		if alpha is not None:
			#only add a pair if their difference is shorter than half of the shorter baseline
			bl1x, bl1y, bl1z = pair[0][1][0]
			bl2x, bl2y, bl2z = pair[1][1][0]
			bench = min([np.hypot(bl1x, bl1y), np.hypot(bl2x, bl2y)])*alpha
			if np.hypot(bl1x-bl2x, bl1y-bl2y) < bench:
				combs.append(pair)
		else:
			combs.append(pair)
	return combs


def get_sub_combs(infile, combs, mode='pm', num=100):
	"""return top num number of combs entries as sorted by 'pm' or 'p'"""
	df = pd.read_csv(infile)
	df['peakmult'] = df['peak']*df['mult']
	df_sorted = pd.DataFrame()
	if mode == 'pm':
		df_sorted = df.sort_values('peakmult', ascending=False)
	elif mode == 'p':
		df_sorted = df.sort_values('peak', ascending=False)
	else:
		raise exception('mode not supported')
	num = min(num, len(combs))
	S = df_sorted.head(n=num).index.tolist()

	return [combs[s] for s in S]

def bl_length(blcoords):
	blx, bly,blz = blcoords
	return np.sqrt(blx*blx + bly*bly + blz*blz)

def run_opp(i, comb, outfile, equiv=None, quiet=False):
	#comb=(('40', (44, 26)), ('41', (44, 38)))
	#i is the index of comb in combs
	file = open(outfile, 'a')
	label1, label2 = comb[0][0], comb[1][0]

	bl1coords, bl2coords = comb[0][1][0], comb[1][1][0]
	bl1L = bl_length(bl1coords); bl2L = bl_length(bl2coords)
	multiplicity = comb[0][1][1]*comb[1][1][1]
	if label1 == label2 and equiv is not None:
		line = ', '.join([str(i), label1,label2,str(0.),str(equiv),str(multiplicity), str(bl1L), str(bl2L)])
		file.write(line+'\n')
		file.close()
		return
	else:
		peak,dT = WS.opp(bl1coords=bl1coords,bl2coords=bl2coords)
		line = ', '.join([str(i), label1,label2,str(dT[0]),str(np.abs(peak[0])),str(multiplicity), str(bl1L), str(bl2L)])
		if not quiet:
			print line
		file.write(line+'\n')
		file.close()


def execute(combsname=None): #for profiling
	CAL = 'psa6622_v003'
	NANTS = 128
	version = 128
	ENTRIES = 1000
	FIRST = 'HERA_{}_all.csv'.format(NANTS)
	SECOND = 'HERA_{0}_pm{1}.csv'.format(NANTS, ENTRIES)
	SECONDm = 'HERA_{0}_p{1}.csv'.format(NANTS, ENTRIES)
	FILE = "../calfiles/HERA_antconfig/antenna_positions_{}.dat".format(version)
	
	#FIRST = 'first.csv'
	
	if combsname:
		print "Loading combs dictionary", combsname
		combs = load_obj(combsname)
	else:
		print "Getting group bl dictionary"
		top_dict, blgps = get_plotsense_dict(file=FILE, NANTS=NANTS)
		print "Looking for appropriate combinations of baselines"
		combs = get_bl_comb(top_dict, alpha=1.)
		save_obj('HERA_{}_combs'.format(NANTS), combs)
	
	if True:
		DT = 0.01
		T1=np.arange(2456681.3, 2456681.7, DT)
		fqs = np.array([.15])
		print 'Starting survey of all baselines'
		global WS 
		WS = w_opp.OppSolver(fqs=fqs, cal=CAL, T1=T1, beam='HERA')
		file = open(FIRST, 'w')
		file.write(',sep,sep2,dT,peak,mult,bl1,bl2\n')
		file.close()

		NJOBS = 12
		print 'Starting Opp with %d instances on %d jobs; dT= %f' % (len(combs), NJOBS, DT)
		Parallel(n_jobs=NJOBS)(delayed(run_opp)(i, comb, FIRST, quiet=True) for i, comb in enumerate(combs))

	if True:
		DT = 0.001
		T1=np.arange(2456681.3, 2456681.7, DT)
		fqs = np.array([.15])
		print '##### search over selected baselines  dT= %f ###'% DT
		global WS
		WS = w_opp.OppSolver(fqs=fqs, cal=CAL, T1=T1, beam='HERA')
		subcombs = get_sub_combs(FIRST, combs, mode='pm', num=ENTRIES)

		
		file = open(SECOND, 'w')
		file.write(',sep,sep2,dT,peak,mult,bl1,bl2\n')
		file.close()

		NJOBS = 4
		print 'Starting Opp with %d instances on %d jobs;' % (len(subcombs), NJOBS)
		#print ',sep,sep2,dT,peak,mult'
		Parallel(n_jobs=NJOBS)(delayed(run_opp)(i, comb, SECOND, equiv=12127.9726, quiet=True) 
			for i, comb in enumerate(subcombs))

	if True:
		DT = 0.001
		T1=np.arange(2456681.3, 2456681.7, DT)
		fqs = np.array([.15])
		print '##### search over selected baselines  dT= %f ###'% DT
		global WS
		WS = w_opp.OppSolver(fqs=fqs, cal=CAL, T1=T1, beam='HERA')
		subcombs = get_sub_combs(FIRST, combs, mode='p', num=ENTRIES)

		
		file = open(SECONDm, 'w')
		file.write(',sep,sep2,dT,peak,mult\n')
		file.close()

		NJOBS = 4
		print 'Starting Opp with %d instances on %d jobs;' % (len(subcombs), NJOBS)
		#print ',sep,sep2,dT,peak,mult'
		Parallel(n_jobs=NJOBS)(delayed(run_opp)(i, comb, SECONDm, equiv=12127.9726, quiet=True) 
			for i, comb in enumerate(subcombs))



	
if __name__=="__main__":
	execute()
	#import IPython; IPython.embed()
