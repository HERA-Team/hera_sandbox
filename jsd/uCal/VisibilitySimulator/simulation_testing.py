import aipy as a
import numpy as np
import optparse, sys, os
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import healpy as hp
import matplotlib.pyplot as plt

#%%
meta,_,mdl,_ = omni.from_npz('./SimulatedVisibilities/simple.2456679.35658.xx.npz')

meta2,_,mdl2,_ = omni.from_npz('../Data/zen.2456679.35658.xx.npz')
#%%
data = mdl['xx']
data2 = mdl2['xx']

bls = [(0, 103), (1, 4), (0, 101), (0, 62), (0, 100), (1, 13), (1, 70), (1, 56), (1, 71), (1, 59), (0, 97), (12, 43), (9, 71), (9, 59), (57, 64)]
conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]

aa = a.cal.get_aa('psa6622_v003', 0.1/203, 0.1, 203)
blLengths = [np.abs((np.round(aa.get_baseline(blIndex[0],blIndex[1]) * a.const.c / 1.5e12))[0]) for blIndex in bls]

freqs = aa.get_afreqs()









plt.figure(1); plt.clf()

#for bl in bls:
#    plt.plot(data[bl][:,0])

bl = bls[12]

plt.subplot(211)
plt.imshow(np.angle(data[bl][:,40:160]))
plt.colorbar()
plt.subplot(212)
plt.imshow(np.angle(data2[bl][:,40:160]))
plt.colorbar()

plt.show()
#%%

plt.figure(2); plt.clf()
test1 = data[bls[3]][:,140]
test2 = data[bls[4]][:,71]
#test1 = data[bls[7]][:,148]
#test2 = data[bls[8]][:,109]
#plt.plot(np.abs(test1)/np.mean(np.abs(test1)))
#plt.plot(np.abs(test2)/np.mean(np.abs(test2)))
plt.plot(np.real(test2/test1))

#%%
print freqs[148]/freqs[109]/(blLengths[8]/blLengths[7])
