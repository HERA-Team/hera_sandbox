#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

MX = -1
DRNG = 2.5
seps = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<1,13> xx',  #  6
    '<1,70> xx',  #  7
    '<1,56> xx',  #  8
    '<1,71> xx',  #  9
    '<1,59> xx',  # 10
    '<0,97> xx',  # 11
    '<12,43> xx', # 12
    '<9,71> xx',  # 13
    '<9,59> xx',  # 14
    '<57,64> xx', # 15
]

conj = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<0,97> xx', # 11
    '<12,43> xx', # 12 
    '<57,64> xx' # 15
]

d,w = {}, {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    for sep in seps:
        if sep in conj: d[sep] = d.get(sep,[]) + [npz[sep].conj()]
        else: d[sep] = d.get(sep,[]) + [npz[sep]]
for sep in seps: d[sep] = n.concatenate(d[sep])
fq1 = n.linspace(.1,.2,d.values()[0].shape[1])
ch1 = n.arange(fq1.size)
sdf = fq1[1]-fq1[0]
if False: # full cross-mult
    d_cat = n.concatenate([d[k] for k in seps], axis=1)
    Cov = n.dot(d_cat.T, d_cat.conj())
    p.subplot(121); C.arp.waterfall(Cov, mx=0, drng=3); p.colorbar()
    p.subplot(122); C.arp.waterfall(Cov, mode='phs'); p.colorbar(); p.show()
constraints = []

for i,s1 in enumerate(seps):
    for j,s2 in enumerate(seps):
        if j <= i: continue
        #if i < 4 or j < 4: continue
        #if i < 10 or j < 10: continue
        print i,s1, j, s2
        cov = n.dot(d[s1].T, d[s2].conj())
        c11 = n.dot(d[s1].T, d[s1].conj())
        c22 = n.dot(d[s2].T, d[s2].conj())
        b1,b2 = float(i+1),float(j+1) # proxy for baseline length
        fq2 = fq1 * b1/b2
        ch2 = n.around((fq2 - .1) / sdf).astype(n.int)
        inds = n.where(n.logical_and(ch2 >= 0, ch2 < fq1.size))[0]
        if len(inds) < 4: continue
        i1,i2 = ch1[inds],ch2[inds]
        #print i1
        #print i2
        #C.arp.waterfall(cov, drng=3); p.show()
        d1 = cov[i1,i2] / c22[i2,i2]
        d2 = n.conj(cov[i1,i2] / c11[i1,i1])
        constraints.append((i1,i2,d1,d2))
        #p.subplot(211); p.plot(n.abs(diag))
        #p.subplot(212); p.plot(n.angle(diag))
#p.show()

CH1 = 30
CH2 = 175
bp1 = n.zeros(fq1.size, dtype=n.complex); bp1[CH1] = 5*1e3 * n.exp(1j*0.8)
bp2 = n.zeros(fq1.size, dtype=n.complex); bp2[CH2] = 1.3 * 1e3
wgt1 = n.zeros(fq1.size); wgt1[CH1] = 1e3
wgt2 = n.zeros(fq1.size); wgt2[CH2] = 1e3
for i in range(4):
    for i1,i2,d1,d2 in constraints:
        for _i1,_i2,_d1,_d2 in zip(i1,i2,d1,d2):
            if wgt2[_i1] > 0 and not n.isnan(_d2) > 0:
                bp2[_i2] += bp2[_i1]/wgt2[_i1] * _d2
                wgt2[_i2] += 1
                print _i1,_i2,_d2
            #elif wgt[_i2] > 0:
            if wgt1[_i2] > 0 and not n.isnan(_d1):
                bp1[_i1] += bp1[_i2]/wgt1[_i2] * _d1
                wgt1[_i1] += 1
                print _i2,_i1, _d1
p.subplot(211); p.plot(n.abs(bp1/wgt1.clip(1,n.inf)))
p.subplot(211); p.plot(n.abs(bp2/wgt2.clip(1,n.inf)))
p.subplot(212); p.plot(n.angle(bp1/wgt1.clip(1,n.inf)))
p.subplot(212); p.plot(n.angle(bp2/wgt2.clip(1,n.inf)))
p.show()
