#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P
import sys

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]

bls = {}
conj = {}
NCOL = len(A_)
for sep in range(1,NCOL):
    for row in [A_,B_,C_,D_]:
        bls[sep] = bls.get(sep,[]) + [(row[i],row[i+sep]) for i in range(NCOL-sep)]
for sep in bls: conj[sep] = [i>j for i,j in bls[sep]]

strbls = {}
conj_bl = {}
for sep in bls:
    strbls[sep] = []
    for (i,j),c in zip(bls[sep],conj[sep]):
        if c: i,j = j,i
        strbls[sep].append('%d_%d' % (i,j))
        conj_bl[ij2bl(i,j)] = c

uv = a.miriad.UV(sys.argv[-1])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

SEP = 1
print ','.join(strbls[SEP])

times, d, f = C.arp.get_dict_of_uv_data(sys.argv[1:], ','.join(strbls[SEP]), 'xx')
for bl in d:
    i,j = bl2ij(bl)
    #print (i,j), conj_bl[bl]
    if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])
    if i == 8 or j == 8: d[bl] = -d[bl]  # XXX remove this line once correct script is run

w = {}
for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)

cal_bl = ij2bl(1,B_[SEP])
bls = [bl for bl in d if bl != cal_bl]
S = 1
for bl in bls:
    taus = []
    for t in range(0,d[cal_bl].shape[0],S):
        print bl2ij(bl),
        g,tau = C.arp.redundant_bl_cal(d[cal_bl][t:t+S],w[cal_bl][t:t+S],d[bl][t:t+S],w[bl][t:t+S],fqs,use_offset=False)
        #P.subplot(311); P.plot(fqs, n.abs(g))
        #P.subplot(312); P.plot(fqs, n.angle(g))
        print 'Delay:', tau
        taus.append(tau)
    P.plot(taus)
P.show()
