import aipy as a, numpy as np, capo as C, pylab as plt
import w_opp
from joblib import Parallel, delayed
import itertools
"""
Computes multiplicities and sensitivities of various baseline combinations
Prints output to OUT
output file should be plotted with scatter_sens.py
"""
ant_dict1 = {'11':(0,26),'12':(0,38),'13':(0,50),'14':(0,23),'15':(0,59),'-11':(0,7),'-12':(0,14)}
ant_dict2 = {'20':(44,26),'21':(44,38),'22':(44,50),'23':(44,23)}
ant_dict3 = {'30':(62,26),'31':(62,38),'32':(62,50)}
ant_dict4 = {'40':(16,26),'41':(16,38)} 
a128_dict1 = {'10':(103,26),'11':(103,38),'12':(103,50),'13':(103,23),'14':(103,59),'1-1':(103,46),'1-2':(95,46)}
a128_dict2 = {'20':(0,26),'21':(0,38),'22':(0,50),'23':(0,23)}
a128_dict3 = {'30':(102,26),'31':(102,38),'32':(102,50)}
a128_dict4 = {'40':(44,26),'41':(44,38)}
#dicts = [ant_dict1,ant_dict2,ant_dict3,ant_dict4]
dicts = [a128_dict1,a128_dict2,a128_dict3,a128_dict4]
CAL = 'psa6622_v003'
OUT = 'corr_res.csv'
file = open(OUT, 'w')
file.write(',sep,sep2,dT,peak,mult\n')
file.close()
def mult(m,n,mm,nn, cal=128):
	if cal==128:
		base = (16-abs(m))*(7-abs(n))*(16-abs(mm))*(7-abs(nn))
	else:
		base = (8-abs(m))*(8-abs(n))*(8-abs(mm))*(8-abs(nn))
	if abs(n) != abs(nn):
		base *= 2
	return base
def get_bl_comb(dicts):
	combs = []
	for dic in dicts:
		for sing in itertools.combinations(dic.iteritems(), 1):
			combs.append((sing[0], sing[0]))
		for pair in itertools.combinations(dic.iteritems(), 2):
			if not(pair[0][1]<0 and pair[1][1]<0):
				combs.append(pair)
	return combs
combs = get_bl_comb(dicts)


WS = w_opp.OppSolver(fqs=np.array([.15]), T1=np.arange(2456681.3, 2456681.7, 0.001), cal=CAL)
print ',sep,sep2,dT,peak,mult'
def run(i, comb):
	#comb=(('40', (44, 26)), ('41', (44, 38)))
	file = open(OUT, 'a')
	key, key2 = comb[0][0], comb[1][0]
	try: m = int(key[0]); n = int(key[1])
	except(ValueError):
		#import IPython; IPython.embed()
		m = int(key[0]); n = int(key[1:])
	try: mm = int(key2[0]); nn = int(key2[1])
	except(ValueError):
		mm = int(key2[0]); nn = int(key2[1:])
	peak,dT = WS.opp(bl1=comb[0][1],bl2=comb[1][1])
	line = ', '.join([str(i), comb[0][0],comb[1][0],str(dT[0]),str(np.abs(peak[0])),str(mult(m,n,mm,nn))])
	print line
	file.write(line+'\n')
	file.close()

Parallel(n_jobs=4)(delayed(run)(i, comb) for i, comb in enumerate(combs))


