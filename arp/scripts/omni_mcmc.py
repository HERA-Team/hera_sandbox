#! /usr/bin/env python

import omnical, aipy, numpy as np, capo
import os, sys, pylab as plt

freqs = np.linspace(.1,.2,1024)

def dtau(g):
    gn = {}
    for pol in g:
        gn[pol] = {}
        for i in g[pol]:
            tau = np.random.normal(scale=.5)
            phs = np.exp(2j*np.pi*tau*freqs)
            phs.shape = (1,-1)
            gn[pol][i] = g[pol][i] * phs
    return gn

def copy(d):
    dn = {}
    for k in d:
        dn[k] = {}
        for j in d[k]: dn[k][j] = d[k][j].copy()
    return dn
            

filename = sys.argv[-1]
fc = '%s.fc.npz' % filename
_,g0,_,_ = capo.omni.from_npz(fc)

aa = aipy.cal.get_aa('hsa7458_v000_HH',np.array([.15]))
ex_ants = [81]
info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants)
reds = info.get_reds()
_,d,f = capo.miriad.read_files([filename], antstr='cross', polstr='xx', decimate=60)


#plt.plot(g0['x'][112].T)
#for i in xrange(10):
for i in xrange(1):
    if i == 0: g = copy(g0)
    else: g = dtau(copy(g0))
    m1,g1,v1 = capo.omni.redcal(d, info, gains=g, removedegen=False)
    m2,g2,v2 = capo.omni.redcal(d, info, gains=g1, vis=v1, uselogcal=False, removedegen=False)
    #plt.plot(g1['x'][112].T)
    for j in g2['x']:
        plt.plot(np.unwrap(np.angle(g0['x'][j][0])))
        plt.plot(np.unwrap(np.angle(g2['x'][j][0])))
plt.show()

import IPython; IPython.embed()
