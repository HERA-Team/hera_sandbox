from w_opp import *
import itertools
fqs = np.linspace(0.145, 0.155, 20)
T1 = np.arange(2456681.45, 2456681.55, 0.0001)

WS = OppSolver(T1=T1, fqs=np.array([.15]))
maxres, T1ac = WS.opp(bl1=(0,26), bl2=(0,38), rephase=0,delay=False, return_series=False)
T1ac += 2456681.5
#Tbox = np.arange(T1ac-0.1, T1ac+0.1, 0.001)
Tbox = T1
WS = OppSolver(T1=Tbox, fqs=fqs)
res = WS.opp(bl1=(0,26), bl2=(0,38), rephase=0,delay=False, return_series=True)
res = np.array(res).T

if False:
	dims = np.ones_like(fqs)
	dims[0] = res.shape[1]
	img = res[0].reshape(dims)
	for i, series in enumerate(res[1:]):
		dims = np.ones_like(fqs)
		dims[i+1] = res.shape[1]
		img = img + series.reshape(dims)
	inds = np.unravel_index(np.argmax(np.abs(img)), img.shape)
	dT = np.array([Tbox[ind]-2456681.5 for ind in inds])
else:
	inds = [np.argmax(series.real) for series in res]
	dT = np.array([Tbox[ind]-2456681.5 for ind in inds])
print dT
import matplotlib.pyplot as plt
plt.figure()
plt.plot(fqs, dT)

import IPython; IPython.embed()
# for pair in itertools.combinations(res, 2):
# 	img = pair[0][:,np.newaxis] + pair[1][np.newaxis,:]
# 	print np.unravel_index(np.argmax(np.abs(img)), img.shape)

# import IPython; IPython.embed()
