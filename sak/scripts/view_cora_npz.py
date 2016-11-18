import numpy as np, healpy as hp
from matplotlib import pylab

def bin2freq(bin): return 100.+(bin/203.)*100.

s = 'all'
spol = ['I','Q','U','V']
freqbin=100
pylab.figure(1)
pylab.clf()

for i,p in enumerate(spol):
    d = np.load(s+p+'.npz')
    hp.mollview(d['map'][freqbin,:],sub=(2,2,i+1),title='%s %f MHz'%(p,bin2freq(freqbin)))
pylab.show()