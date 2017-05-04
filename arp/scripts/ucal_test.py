#! /usr/bin/env python
import numpy as np, pylab as plt, aipy, capo, omnical
import capo.hex

aa = aipy.cal.get_aa('hsa7458_v000_HH', np.array([.15]))
capo.hex.aa_to_info_hera(aa)
info = capo.hex.aa_to_info_hera(aa)
reds = info.get_reds()

POL = 'yy'
files = ['zen.2457678.39660.%s.HH.uvc' % POL]
blgp1 = [((9,43),(88,43)), ((53,97),(80,97)), ((72,31),(72,96))]
blgp2 = [((20,43),(9,43)),((20,97),(53,97)),((72,20),(72,31))]
bls = blgp1 + blgp2
antstr = ','.join(['%d_%d,%d_%d' % (bl[0] + bl[1]) for bl in bls])

print 'Reading data for', files
fm,fg,fv,fx = capo.omni.from_npz(files[0]+'.fc.npz')
m,g,v,x = capo.omni.from_npz(files[0][:-len('HH.uvc')]+'npz')
info,data,flgs = capo.miriad.read_files(files, antstr=antstr, polstr=POL)

ALL_ANTS = g[POL[0]].keys()
NTIMES,NCHAN = g[POL[0]][ALL_ANTS[0]].shape
ANTS = [88,9,20,89,43]
rfi = capo.xrfi.omni_chisq_to_flags(m['chisq']) # slow
#import sys; sys.exit(0)
#import IPython; IPython.embed()

print 'Solving for degeneracies'
phsgrad = np.empty((NTIMES,3,NCHAN), np.float)
for t in xrange(NTIMES):
    gdata = np.angle([g[POL[0]][i][t]/fg[POL[0]][i].flatten() for i in ALL_ANTS])
    d = {}
    for i,ai in enumerate(ALL_ANTS):
        d['%f*dphsdx+%f*dphsdy+%f*dphsdz' % tuple(aa[ai].pos)] = gdata[i]
    ls = capo.linsolve.LinearSolver(d)
    sols = ls.solve()
    phsgrad[t] = np.array([sols['dphsdx'], sols['dphsdy'], sols['dphsdz']])

print 'Removing degeneracies from omnical solutions'
pg = {POL[0]:{}}
apg = {POL[0]:{}}
for i in ALL_ANTS: pg[POL[0]][i] = g[POL[0]][i] / np.exp(1j*np.dot(aa[i].pos, phsgrad))
avgphs = np.average([np.angle(pg[POL[0]][i]/fg[POL[0]][i].flatten()) for i in ALL_ANTS], axis=0)
avgamp = np.average([np.abs(pg[POL[0]][i]/fg[POL[0]][i].flatten()) for i in ALL_ANTS], axis=0)
for i in ALL_ANTS: apg[POL[0]][i] = pg[POL[0]][i] / np.exp(1j*avgphs) / avgamp

apv = {POL:{}}
for bl in v[POL]:
    ai,aj = bl
    dphs = np.dot(aa[ai].pos-aa[aj].pos, phsgrad)
    apv[POL][bl] = v[POL][bl] * np.exp(1j*dphs) * avgamp**2

for bl in data:
    ai,aj = bl
    data[bl][POL] /= apg[POL[0]][ai] * apg[POL[1]][aj].conj()

# begin ucal stuff
NCHAN = 256
fq = np.linspace(.1,.2,NCHAN)
sdf = fq[1] - fq[0]
blpairs = blgp1 #+ blgp2
sols = {}
chpairs = {}
eqs = {}
for bl1,bl2 in blpairs:
    d1 = aa[bl1[0]].pos - aa[bl1[1]].pos
    d2 = aa[bl2[0]].pos - aa[bl2[1]].pos
    r = np.around(np.sqrt(np.dot(d1,d1)) / np.sqrt(np.dot(d2,d2)), 3)
    #r = .75
    try: apv1 = data[bl1][POL]
    except(KeyError): apv1 = data[bl1[::-1]][POL].conj()
    try: apv2 = data[bl2][POL]
    except(KeyError): apv2 = data[bl2[::-1]][POL].conj()
    apv1.shape = (apv1.shape[0],NCHAN,-1); apv1 = np.average(apv1, axis=-1)
    apv2.shape = (apv2.shape[0],NCHAN,-1); apv2 = np.average(apv2, axis=-1)
    chs = np.array([np.around((f1*r - fq[0])/sdf).astype(np.int) for f1 in fq[:-52/4]])
    valid = np.where(chs > 52 / 4, 1, 0)
    ch1s = np.arange(NCHAN).compress(valid)
    ch2s = chs.compress(valid)
    for ch1,ch2 in zip(ch1s,ch2s):
        ublkey = '0'.join(map(lambda x: str(int(x)), np.around(fq[ch1]*d1*10, 0)))
        ublkey = ublkey.replace('-','n')
        eqs['bp%d * ubl%s' % (ch1,ublkey)] = apv1[:2,ch1]
        eqs['bp%d * ubl%s' % (ch2,ublkey)] = apv2[:2,ch2]
    sols[bl1] = apv1[:,ch1s] / apv2[:,ch2s]
    chpairs[bl1] = (ch1s, ch2s)
#plt.ylim(-.5,.5)

#wgts = {}
#for k in eqs: wgts[k] = np.ones(eqs[k].shape, dtype=np.float)
print len(eqs)
ls = capo.linsolve.LogProductSolver(eqs)

A = ls.ls_phs.get_A()
AtA = A[...,0].T.dot(A[...,0])
U,S,V = np.linalg.svd(AtA)
plt.plot(S); plt.show()
bp_prms = [k for k in ls.ls_phs.prm_order.keys() if k.startswith('bp')]
def mycmp(x,y): return cmp(int(x[2:]),int(y[2:]))
bp_prms.sort(mycmp)
bp_ind = np.array([ls.ls_phs.prm_order[k] for k in bp_prms])
capo.plot.waterfall(U[:,bp_ind], mode='lin'); plt.show()
plt.plot(U[:10,bp_ind].T); plt.show()

import IPython; IPython.embed()

sols = {}
bps = {}
print 'Solving'
sols[0] = ls.solve()
print 'Finished log'
def mycmp(x,y): return cmp(int(x[2:]),int(y[2:]))
ch = [k for k in sols[0].keys() if k.startswith('bp')]
ch.sort(mycmp)
bps[0] = np.array([sols[0][k] for k in ch]).T
for i in xrange(1,7):
    ls = capo.linsolve.LinProductSolver(eqs, sols[i-1])
    print 'Solving', i
    sols[i] = ls.solve()
    print 'Finished', i
    bps[i] = np.array([sols[i][k] for k in ch]).T
import IPython; IPython.embed()
