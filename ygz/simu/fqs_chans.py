from w_opp import OppSolver
import aipy as a, numpy as np, capo as C, pylab as plt, seaborn as sns
"""
This script uses the w_opp module to plot the effect of number of frequency channels 
"""

Afqs = np.linspace(.1, .2, 203)
FQL = [np.array([.15]), Afqs[125:146], Afqs[110:161]]
bp = np.load('../64code/bandpass.npz')['bandpass']
T1 = np.arange(2456681.4, 2456681.60, 0.002)
RES = np.zeros((len(FQL), T1.size), dtype=np.complex)
for i, fqs in enumerate(FQL):
	WS = OppSolver(fqs=fqs, bandpass=bp, T1=T1)
	res = WS.opp(bl1=(0,103), bl2=(0,103), return_series=False)[0]
	res2 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True)
	RES[i,:] = res2/np.abs(res)

plt.figure()
for i, fqs in enumerate(FQL):
	plt.plot(T1-WS.T0, np.abs(RES[i]), label=str(fqs.size))
plt.legend()
plt.show()

import IPython; IPython.embed()