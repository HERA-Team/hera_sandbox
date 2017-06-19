from w_opp import *
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
N = 100
T1 = np.arange(2456681.4, 2456681.6, 0.001)
WS = OppSolver(fqs=np.linspace(.14,.16,10))
S_simu = Parallel(n_jobs=8)(delayed(WS.opp)(bl1=(0,26), bl2=(0,26), sky=True, return_series=True) for i in xrange(N))

S_simu = np.means(np.asarray(S_simu), axis=0)

S_anal = WS.opp(bl1=(0,26), bl2=(0,26), sky=False, return_series=True)


plt.figure()
plt.plot(T1, S_simu, label='simu')
plt.plot(T1, S_anal, label='anal')
plt.show()