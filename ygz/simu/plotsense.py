import aipy as a, numpy as np, capo as C, pylab as plt
import w_opp

ant_dict1 = {'11':(0,26),'12':(0,38),'13':(0,50),'14':(0,23),'15':(0,59),'-11':(0,7),'-12':(0,14)}
ant_dict2 = {'20':(44,26),'21':(44,38),'22':(44,50),'23':(44,23)}
ant_dict3 = {'30':(62,26),'31':(62,38),'32':(62,50)}
ant_dict4 = {'40':(16,26),'41':(16,38)} 
dicts = [ant_dict1,ant_dict2,ant_dict3,ant_dict4]

def mult(m,n,mm,nn):
	return (8-abs(m))*(8-abs(n))*(8-abs(mm))*(8-abs(nn))
print 'sep sep2 dT peak mult'
finished = []
for dic in dicts:
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
			if key == key2: 
				print key,key2,0.,1.,mult(m,n,mm,nn)
			else:
				peak,dT = w_opp.w_opp(dic[key],dic[key2])
				print key,key2,dT,np.abs(peak),mult(m,n,mm,nn)
				finished.append((key,key2))


