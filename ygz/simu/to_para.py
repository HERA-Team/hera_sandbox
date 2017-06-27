from w_opp import *
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
RPHS = 'auto'
N = 16
T1 = np.arange(2456681.4, 2456681.6, 0.001)
global WS
WS = OppSolver(fqs=np.linspace(.14,.16,10))

def run_simu(bl1=(0,26), bl2=(0,26)):
	return WS.opp(bl1=bl1, bl2=bl2, sky=True, return_series=True, rephase=RPHS)
S_simu = Parallel(n_jobs=16)(delayed(run_simu)(bl1=(0,26), bl2=(0,38)) for i in xrange(N))
#S_simu = [WS.opp(bl1=(0,26), bl2=(0,26), sky=True, return_series=True) for i in xrange(N)]
S_simu = np.mean(np.asarray(S_simu), axis=0)

S_anal = WS.opp(bl1=(0,26), bl2=(0,38), sky=False, return_series=True, rephase=RPHS)
print np.abs(S_simu)/np.abs(S_anal)
S_anal /= np.amax(np.abs(S_anal))
S_simu /= np.amax(np.abs(S_simu))
#import IPython; IPython.embed()
plt.figure()
plt.plot(T1, S_simu.real, label='simu')
plt.plot(T1, S_anal.real, label='anal')
plt.legend()
plt.show()
import IPython; IPython.embed()