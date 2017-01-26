import aipy as a, numpy as np, capo as C, pylab as plt
import w_opp
from joblib import Parallel, delayed


ant_dict1 = {'11':(0,26),'12':(0,38),'13':(0,50),'14':(0,23),'15':(0,59),'-11':(0,7),'-12':(0,14)}
ant_dict2 = {'20':(44,26),'21':(44,38),'22':(44,50),'23':(44,23)}
ant_dict3 = {'30':(62,26),'31':(62,38),'32':(62,50)}
ant_dict4 = {'40':(16,26),'41':(16,38)} 
a128_dict1 = {'10':(103,26),'11':(103,38),'12':(103,50),'13':(103,23),'14':(103,59),'-11':(103,12),'-12':(103,54)}
a128_dict2 = {'20':(0,26),'21':(0,38),'22':(0,50),'23':(0,23)}
a128_dict3 = {'30':(102,26),'31':(102,38),'32':(102,50)}
a128_dict4 = {'40':(44,26),'41':(44,38)}
#dicts = [ant_dict1,ant_dict2,ant_dict3,ant_dict4]
dicts = [a128_dict1,a128_dict2,a128_dict3,a128_dict4]
CAL = 'psa6622_v003'
OUT = 'corr_res.csv'
file = open(OUT, 'w')
file.write('sep sep2 dT peak mult\n')
file.close()
def mult(m,n,mm,nn, cal=128):
	if cal==128:
		return (16-abs(m))*(7-abs(n))*(16-abs(mm))*(7-abs(nn))
	else:
		return (8-abs(m))*(8-abs(n))*(8-abs(mm))*(8-abs(nn))

finished = []
WS = w_opp.OppSolver(fq=.15, cal=CAL)
print 'sep sep2 dT peak mult'
#def run(key, key2, dic):
def run(dic):
	file = open(OUT, 'a')
	for key in dic.keys():
		for key2 in dic.keys():
			if (key,key2) in finished: continue
			try: m = int(key[0]); n = int(key[1])
			except(ValueError):
				#import IPython; IPython.embed()
				m = int(key[:2]); n = int(key[2])
			try: mm = int(key2[0]); nn = int(key2[1])
			except(ValueError):
				mm = int(key2[:2]); nn = int(key2[2])
			# if key == key2: 
			# 	print key,key2,0.,1.,mult(m,n,mm,nn)
			# else:
			peak,dT = WS.w_opp(dic[key],dic[key2])
			line = ' '.join([str(key),str(key2),str(dT),str(np.abs(peak)),str(mult(m,n,mm,nn))])
			print line
			file.write(line+'\n')
			finished.append((key,key2))
	file.close()
Parallel(n_jobs=len(dicts))(delayed(run)(dic) for dic in dicts)


