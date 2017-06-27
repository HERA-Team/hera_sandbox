#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob

SEPS = [
    (0,103) , #  1
    #(1,4) ,   #  2
    (0,26), # 2 for psa128 v3
    (0,101) , #  3
    (0,62) ,  #  4
]
'''
    (0,100) , #  5
    (1,13) ,  #  6
    (1,70) ,  #  7
    (1,56) ,  #  8
    (1,71) ,  #  9
    (1,59) ,  # 10
    (0,97) ,  # 11
    (12,43) , # 12
    (9,71) ,  # 13
    (9,59) ,  # 14
    (57,64) , # 15
]'''

#fqs = np.linspace(.1,.2,bandpass.size)

sets = {
    #'day0' : sys.argv[1:],
    'day0' : glob.glob('zen.2456680.*.xx.npz'),
    'day1' : glob.glob('zen.2456681.*.xx.npz'),
    'day2' : glob.glob('zen.2456682.*.xx.npz'),
    #'day3' : glob.glob('zen.2456683.*.xx.npz'),
    'day4' : glob.glob('zen.2456684.*.xx.npz'),
    'day5' : glob.glob('zen.2456685.*.xx.npz'),
}

data,wgts = {}, {}
lsts, chisqs = {}, {}
for s in sets:
    if not lsts.has_key(s):
        meta, gains, vismdl, xtalk = capo.omni.from_npz(sets[s], bls=SEPS, verbose=True)
        lsts[s] = meta['lsts']
    chisqs[s] = meta['chisq']
    for pol in vismdl:
        for bl in vismdl[pol]:
            k = (s,pol,bl)
            data[k] = vismdl[pol][bl]
            #if bl in CONJ: data[k] = data[k].conj()
            #data[k] *= bandpass[:,CH0:CH0+NCHAN]
            wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)
            #wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1./chisqs[s])

lsts_g = {}
data_g,wgts_g = {}, {}
if False: # new expensive gridding
    for k in data:
        print 'Gridding', k
        lsts_g[k],data_g[k],wgts_g[k] = capo.oqe.lst_grid(lsts[k[0]], data[k], wgts=wgts[k])
else: # old cheap gridding (nearest)
    for k in data:
        print 'Gridding', k
        lsts_g[k],data_g[k],wgts_g[k] = capo.oqe.lst_grid_cheap(lsts[k[0]], data[k], wgts=wgts[k], lstbins=3000)
        #lst_res = np.average(lsts[k[0]][1:] - lsts[k[0]][:-1])
        #inds = capo.oqe.lst_align(lsts, lstres=lst_res)
        #data_g[k],wgts_g[k],lsts_g[k] = capo.oqe.lst_align_data(inds, dsets=data, wgts=wgts, lsts=lsts[k[0]])
        #for s in sets: chisqs[s] = chisqs[s][inds[s]].T
    

def k2eq(k): return 'g%s * bl%d_%d' % ((k[0],) + k[-1])

data_eqs, wgts_eqs = {}, {}
sol0 = {}
CH = 110
for k in data_g:
    if not k[-1] == (0,103): continue
    data_eqs[k2eq(k)], wgts_eqs[k2eq(k)] = data_g[k][:,CH].copy(), wgts_g[k][:,CH].copy()
    sol0['g'+k[0]] = np.ones_like(data_g[k][:,CH])
    sol0['bl%d_%d' % k[-1]] = data_g[k][:,CH].copy()

ls = capo.linsolve.LinProductSolver(data_eqs, wgts_eqs, sol0)
sol1 = ls.solve()
ls = capo.linsolve.LinProductSolver(data_eqs, wgts_eqs, sol1)
sol2 = ls.solve()

for g in [k for k in sol2.keys() if k.startswith('g')]:
    plt.plot(sol2[g])

plt.show()

import IPython; IPython.embed()
