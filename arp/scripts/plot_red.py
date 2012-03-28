#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P
import sys, scipy

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

USE_CAL = True
IMPROVE = False

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
if USE_CAL: aa = a.cal.get_aa('psa898_v002', fqs)

#seps = ['0,1']
#seps = ['0,2']
seps = ['0,6']
plot_bls = bls[seps[0]]
#plot_bls = [ij2bl(4,15)]
#seps = bls.keys()
strbls = ','.join([strbls[sep] for sep in seps])
print strbls

times, d, f = C.arp.get_dict_of_uv_data(sys.argv[1:], strbls, 'xx', verbose=True)
for bl in d:
    i,j = bl2ij(bl)
    if USE_CAL:
        d[bl] = aa.phs2src(d[bl], 'z', i,j)
        d[bl] /= n.median(aa.passband(i,j))
    else:
        if i == 8 or j == 8: d[bl] = -d[bl]  # XXX remove this line once correct script is run
    if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])

w = {}
for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)

NANT = ANTPOS.size
ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))

d_sum, d_wgt = {}, {}

#calrow,calcol = 0,0
calrow,calcol = 1,0
calant = ANTPOS[calrow,calcol]
cal_bl = {}
for sep in seps: cal_bl[sep] = [bl for bl in bls[sep] if ANTPOS[calrow,calcol] in bl2ij(bl)][0]

def avg(dsum,wsum): return dsum / n.where(wsum == 0, 1, wsum)

for sep in seps:
    cbl = cal_bl[sep]
    i0,j0 = bl2ij(cbl)
    if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
    d_sum[sep] = d[cbl].copy(); d_wgt[sep] = w[cbl].copy()
    for bl in bls[sep]:
        if bl == cbl or not bl in plot_bls: continue
        i,j = bl2ij(bl)
        if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
        if IMPROVE:
            g,tau,info = C.arp.redundant_bl_cal(d[cbl],w[cbl],d[bl],w[bl],fqs,use_offset=False)
            gain = n.median(n.abs(g))
        else:
            g,tau,info = C.arp.redundant_bl_cal(d[cbl],w[cbl],d[bl],w[bl],fqs,use_offset=False,maxiter=0)
            gain = 1
        d_sum[sep] += d[bl] * n.exp(-2j*n.pi*fqs*tau) / gain
        d_wgt[sep] += w[bl]
        P.subplot(211); P.semilogy(fqs, n.abs(g)/gain, label='%d,%d'%(i,j))
        P.subplot(212); P.plot(fqs, n.angle(g), label='%d,%d'%(i,j))
        print (i,j)
        if USE_CAL: print '  Prev:', aa[j].get_params(['dly'])['dly'] - aa[i].get_params(['dly'])['dly']
        print '   Dly:', tau, 'Dtau:', info['dtau']
        print '  Gain:', gain
    P.subplot(211); P.semilogy(fqs, avg(n.abs(d[cbl]).sum(axis=0), w[cbl].sum(axis=0)))
        
P.legend()
P.show()

for sep in seps:
    cbl = cal_bl[sep]
    dcal = avg(d[cbl],w[cbl])
    P.subplot(231); C.arp.waterfall(dcal, drng=2)
    P.subplot(232); C.arp.waterfall(dcal, mode='phs')
    P.subplot(233); C.arp.waterfall(n.fft.fftshift(C.arp.clean_transform(d[cbl],w[cbl]), axes=-1), drng=2.5)
    d = avg(d_sum[sep], d_wgt[sep])
    P.subplot(234); C.arp.waterfall(d, drng=2)
    P.subplot(235); C.arp.waterfall(d, mode='phs')
    P.subplot(236); C.arp.waterfall(n.fft.fftshift(C.arp.clean_transform(d_sum[sep],d_wgt[sep]), axes=-1), drng=2.5)
    P.show()

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
