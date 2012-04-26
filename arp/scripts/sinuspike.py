#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, sys

uv = a.miriad.UV(sys.argv[-1])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)
times, d, f = C.arp.get_dict_of_uv_data(sys.argv[1:], '0_16,1_17,2_18,3_19', 'xx', verbose=True)

def kurtosis(d,w,axis=0,rm_avg=True):
    '''Only works for wgts consisting of 1,0'''
    wgt = n.sum(w, axis=axis).clip(1,n.Inf)
    if rm_avg:
        # XXX this step may fail if axis!=0
        avg = n.sum(d, axis=axis) / wgt
        res = d - avg
    else: res = d
    m4 = n.sum(n.abs(res)**4, axis=axis) / wgt
    m2 = n.sum(n.abs(res)**2, axis=axis) / wgt
    return m4/n.where(m2==0, 1, m2**2) - 3

for bl in d:
    print a.miriad.bl2ij(bl)
    flg,mdl = [], []
    f[bl][:,:140] = 1
    f[bl][:,1920:] = 1
    wgt = n.logical_not(f[bl])
    d[bl] *= wgt
    k = kurtosis(d[bl],wgt)
    p.subplot(211); p.plot(k)
    f[bl][:,n.where(k>3)] = 1
    for i in range(len(times)):
        mi,fi = C.arp.sinuspike(d[bl][i], fqs, f=f[bl][i], nsig=2)
        mi *= 0
        mdl.append(mi)
        flg.append(fi)
    flg = n.array(flg)
    mdl = n.array(mdl)
    wgt = n.logical_not(flg).astype(n.int)
    res = (d[bl] - mdl) * wgt
    kurt = kurtosis(res,wgt,rm_avg=False)
    #p.subplot(211); p.plot(kurt); p.plot(kurtosis(d[bl], n.ones_like(d[bl])))
    p.subplot(211); p.plot(kurt)
    p.subplot(211); p.plot(kurtosis(res,wgt,rm_avg=False,axis=1))
    flg[:,n.where(kurt>3)] = 1
    #p.subplot(131); C.arp.waterfall(d[bl], mx=0, drng=2)
    #p.subplot(132); C.arp.waterfall(  mdl, mx=0, drng=2)
    #p.subplot(133); C.arp.waterfall(  res*n.logical_not(flg), drng=2)
    p.subplot(212); C.arp.waterfall(  res*n.logical_not(flg), drng=2)
    p.show()
