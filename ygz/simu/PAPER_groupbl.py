import numpy as np, pandas as pd
import itertools, pickle
import w_opp, aipy as a
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


def get_plotsense_dict(cal, NTOP=None, NANTS=112):
	aa = a.cal.get_aa(cal, np.array([.15]))
	antpos = [aa.get_baseline(0,i,src='z') for i in xrange(NANTS)]
	antpos = np.array(antpos) * a.const.len_ns / 100.
	#antpos = np.loadtxt(file)
	if NANTS is None:
		NANTS = antpos.shape[0]
	M,N = aa.ant_layout.shape
	bl_gps = {}
	for m in xrange(M):
		for n in xrange(N):
			for mm in xrange(M):
				for nn in xrange(N):
					if m == mm and n == nn: continue
					i, j = aa.ant_layout[m,n], aa.ant_layout[mm,nn]
					#a1pos = antpos[i]
					#a2pos = antpos[j]
					blm, bln = mm-n, nn-n
					if blm <= 0: continue
					#blx, bly, blz = a2pos-a1pos
					has_entry = False
					for key in bl_gps.keys():
						if ( key[0]==blm and key[1]==bln ):
							has_entry = True
							bl_gps[key] = bl_gps.get(key, []) + [(i,j)]
							break
					if not has_entry:
						bl_gps[(blm, bln)] = [(i,j)]

	n_unique = len(bl_gps.keys())
	n_total = NANTS*(NANTS-1)/2
	
	top_dict = {}
	for i, key in enumerate(bl_gps.keys()):
		mult = len(bl_gps[key])
		label = str(key[0])+'_'+str(key[1]) #label = blm, bln
		i,j = bl_gps[key][0]
		a1pos = antpos[i]; a2pos = antpos[j]
		blx, bly, blz = a2pos-a1pos
		top_dict[label] = ((blx, bly, blz), mult) #dict[example_bl] = ((blx, bly, blz), multiplicity)
	return top_dict, bl_gps

def get_bl_comb(top_dict, alpha=None):
	combs = []
	for sing in itertools.combinations(top_dict.iteritems(), 1):
		combs.append((sing[0], sing[0]))
	for pair in itertools.combinations(top_dict.iteritems(), 2):
		# !!!!!!!!!
	
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
	NANTS = 112
	version = 128
	ENTRIES = 1000
	FIRST = 'PAPER_{}_all.csv'.format(version)
	SECOND = 'PAPER_{0}_pm.csv'.format(version)
	SECONDm = 'PAPER_{0}_p.csv'.format(version)
	
	#FIRST = 'first.csv'
	
	if combsname:
		print "Loading combs dictionary", combsname
		combs = load_obj(combsname)
	else:
		print "Getting group bl dictionary"
		top_dict, blgps = get_plotsense_dict(cal=CAL, NANTS=NANTS)
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

		NJOBS = 4
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
