#! /usr/bin/env python

import aipy, numpy as np, pylab as plt, capo
import sys

POL = sys.argv[-1].split('.')[-3]
info, data, flgs = capo.miriad.read_files(sys.argv[1:], antstr='all', polstr=POL)
ants = {}
for i,j in data.keys():
    ants[i] = ants.get(i,len(ants))
    ants[j] = ants.get(j,len(ants))
print len(ants), ants
N = len(ants)
for cnt,ch in enumerate(xrange(100,900,50)):
    M = np.zeros((N,N), dtype=np.complex)
    for i,j in data.keys():
        #if i == j: continue
        x,y = ants[i], ants[j]
        M[x,y] = data[(i,j)][POL][0,ch]
        if i != j: M[y,x] = data[(i,j)][POL][0,ch].conjugate()
    U,S,V = np.linalg.svd(M)
    plt.subplot(4,4,cnt+1); capo.plot.waterfall(U, drng=3)
plt.show()
import IPython; IPython.embed()
