#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys
from arp_paper import *

d, w = {}, {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    f = n.load(filename)
    for bl in f.files:
        i = a.miriad.bl2ij(int(bl))[0]
        #if i <= 16: continue
        _ant = f[bl]
        #_ant = _ant[:,72:910]
        _wgt = n.where(_ant == 0, 0., 1.)
        #sig = n.std(_ant.flatten().compress(_wgt.flatten()))
        #_wgt = n.where(n.abs(_ant) > 2*sig, 0, _wgt)
        #_ant *= _wgt
        d[i] = d.get(i,[]) + [_ant.sum(axis=0)]
        w[i] = w.get(i,[]) + [_wgt.sum(axis=0)]

for i in d:
    d[i] = n.array(d[i])
    w[i] = n.array(w[i])

#_minavg = None
_mindly = None
for i in d:
    p.subplot(221)
    avg = d[i].sum(axis=0) / w[i].sum(axis=0).clip(1,n.Inf)
    p.plot(avg)
    #if _minavg is None: _minavg = n.where(aavg > 0, avg, n.Inf)
    #else: _minavg = n.where(n.logical_and(aavg < n.abs(_minavg), aavg > 0), avg, _minavg)

    p.subplot(222)
    dly = n.fft.fft(avg)
    adly = n.abs(dly)
    p.semilogy(adly)
    if _mindly is None: _mindly = n.where(adly > 0, dly, n.Inf)
    else: _mindly = n.where(n.logical_and(adly < n.abs(_mindly), adly > 0), dly, _mindly)

_sumdly, _sumwgt = 0, 0
for i in d:
    avg = d[i].sum(axis=0) / w[i].sum(axis=0).clip(1,n.Inf)
    dly = n.fft.fft(avg)
    adly = n.abs(dly)
    #_wgt = n.where(adly < 2*n.abs(_mindly), 1., 0.)
    #_wgt = n.where(adly < 4*n.abs(_mindly), 1., 0.)
    _wgt = n.where(adly < 8*n.abs(_mindly), 1., 0.)
    _sumdly += dly * _wgt
    _sumwgt += _wgt
print _sumwgt.max(), n.average(_sumwgt)

def search_steps(spec):
    step = []
    wgt = n.where(n.abs(spec) == 0, 0., 1.)
    for thresh in range(1,spec.size-1):
        left = n.sqrt(n.sum(n.abs(spec[:thresh])**2) / max(1, n.sum(wgt[:thresh])))
        right = n.sqrt(n.sum(n.abs(spec[thresh:])**2) / max(1, n.sum(wgt[thresh:])))
        step.append(left - right)
    return n.array([0] + step + [0])

_sumdly = _sumdly / _sumwgt.clip(1,n.Inf)

p.subplot(2,2,3)
#p.plot(_minavg)
#_sumdly[:131] = 0
#_sumdly[-130:] = 0
#_mindly[:131] = 0
#_mindly[-130:] = 0
_minavg = n.fft.ifft(_mindly)
_sumavg = n.fft.ifft(_sumdly)
#_minavg[:70] = 0
#_minavg[910:] = 0

#p.plot(n.convolve(n.abs(_minavg), n.ones(8), mode='same'))
#p.plot(n.convolve(n.abs(_sumavg), n.ones(8), mode='same'))
#p.plot(n.abs(_minavg))
#p.plot(n.abs(_sumavg))
p.plot(_minavg)
p.plot(_sumavg)
p.plot(search_steps(_minavg))
p.plot(search_steps(_sumavg))

p.subplot(2,2,4)
#p.plot(n.abs(n.fft.fft(_minavg)))
p.plot(n.abs(_mindly))
p.plot(n.abs(_sumdly))

#p.plot(_min)

#tot,totwgt = 0,0
#for bl in ant:
#    tot += ant[bl]
#    totwgt += wgt[bl]


