#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P
import sys

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

bls = {}
conj = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ri,ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                sep = '%d,%d' % (rj-ri, cj-ci)
                bls[sep] = bls.get(sep,[]) + [(ANTPOS[ri,ci],ANTPOS[rj,cj])]
for sep in bls.keys():
    if sep == '0,0' or len(bls[sep]) < 2: del(bls[sep])
for sep in bls:
    conj[sep] = [i>j for i,j in bls[sep]]

strbls = {}
conj_bl = {}
for sep in bls:
    strbls[sep] = []
    bl_list = []
    for (i,j),c in zip(bls[sep],conj[sep]):
        if c: i,j = j,i
        bl_list.append(ij2bl(i,j))
        strbls[sep].append('%d_%d' % (i,j))
        conj_bl[ij2bl(i,j)] = c
    bls[sep] = bl_list
    strbls[sep] = ','.join(strbls[sep])
    print sep, strbls[sep]

uv = a.miriad.UV(sys.argv[-1])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

#seps = ['0,1','1,1','1,0','0,2']
seps = bls.keys()
strbls = ','.join([strbls[sep] for sep in seps])
print strbls

times, d, f = C.arp.get_dict_of_uv_data(sys.argv[1:], strbls, 'xx', verbose=True)
for bl in d:
    i,j = bl2ij(bl)
    if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])
    if i == 8 or j == 8: d[bl] = -d[bl]  # XXX remove this line once correct script is run

w = {}
for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)

NANT = ANTPOS.size
ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))

calrow,calcol = 0,0
calant = ANTPOS[calrow,calcol]
cal_bl = {}
for sep in seps: cal_bl[sep] = [bl for bl in bls[sep] if bl2ij(bl)[0] == ANTPOS[calrow,calcol]][0]

P = n.zeros((3,NANT), dtype=n.float)
M = n.zeros((3,1), dtype=n.float)
P[0,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M[0] = 0
P[1,ANTIND[ANTPOS[calrow,calcol+1]]] = 1e6; M[0] = 0
P[2,ANTIND[ANTPOS[calrow+1,calcol]]] = 1e6; M[0] = 0
print P
    
for sep in seps:
    cbl = cal_bl[sep]
    i0,j0 = bl2ij(cbl)
    if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
    for bl in bls[sep]:
        if bl == cbl: continue
        i,j = bl2ij(bl)
        if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
        Pline = n.zeros((1,NANT), dtype=n.float)
        Pline[0,ANTIND[j]] += 1; Pline[0,ANTIND[i]] += -1; Pline[0,ANTIND[j0]] += -1; Pline[0,ANTIND[i0]] += 1
        Mline = n.zeros((1,1), dtype=n.float)
        g,tau = C.arp.redundant_bl_cal(d[cbl],w[cbl],d[bl],w[bl],fqs,use_offset=False)
        Mline[0,0] = tau
        P = n.append(P, Pline, axis=0)
        M = n.append(M, Mline, axis=0)
        print ''.join(['v-0+^'[int(c)+2] for c in Pline.flatten()]), Mline, (i,j), (i0,j0)

C = n.linalg.lstsq(P,M)[0]
C.shape = ANTPOS.shape
print n.around(C,2)
    #if True:
    #    print bl2ij(bl), 'Delay:', tau
    #    #P.subplot(211); P.plot(fqs, n.abs(g))
    #    #P.subplot(212); P.plot(fqs, n.angle(g))
    #    #P.plot(fqs, n.angle(n.sum(d[cal_bl]*n.conj(d[bl]), axis=0)))
    #else:
    #    taus = []
    #    tau = 0.
    #    for t in range(0,d[cal_bl].shape[0],S):
    #        g,tau = C.arp.redundant_bl_cal(d[cal_bl][t:t+S],w[cal_bl][t:t+S],d[bl][t:t+S],w[bl][t:t+S],fqs,use_offset=False, tau=tau)
    #        print 'Delay:', tau
    #        taus.append(tau)
    #    P.plot(taus)
    
#P.show()
