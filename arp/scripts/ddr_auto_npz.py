#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys
from arp_paper import *

autos = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    autos[filename] = {}
    f = n.load(filename)
    for bl in f.files: autos[filename][int(bl)] = f[bl]

keys = autos.keys(); keys.sort()

ant,wgt = {},{}
__sum, __wgt = {}, {}

_min = None
for cnt, filename in enumerate(keys):
    print filename
    for bl in autos[filename]:
        i = a.miriad.bl2ij(bl)[0]
        #if i != 4: continue
        _ant = autos[filename][bl]
        _ant[:20] = 0; _ant[-20:] = 0
        #_ant[:,:31] = 0; _ant[:,200:] = 0
        _ant = _ant[:,31:200]
        _wgt = n.where(_ant == 0, 0., 1.)
        sig = n.std(_ant.flatten().compress(_wgt.flatten()))
        _wgt = n.where(n.abs(_ant) > 2*sig, 0, _wgt)
        _ant *= _wgt
        _tant = 1e-3 * n.ones_like(_ant); _tant[:,90:] = 0
        #_tant = 1e-3 * n.random.normal(size=_ant.shape); _tant[:,90:] = 0
        #_tant = _ant.copy(); _tant[:,90:] = 0
        _tant *= _wgt
        _dly = n.fft.fft(_ant, axis=1)
        _tdly = n.fft.fft(_tant, axis=1)
        #_ker = n.fft.fft(_wgt, axis=1)
        #g = n.sqrt(n.average(_wgt**2, axis=1))
        #for t in range(_dly.shape[0]):
        #    print _dly[t].sum(), _ker[t].sum(), g[t]
        #    if g[t] == 0: continue
        #    d,info = a.deconv.clean(_dly[t].copy(), _ker[t].copy())
        #    _dly[t] = d + info['res'] / g[t]
        _dly[:,:31] = 0; _dly[:,-30:] = 0
        _dly[:,77:93] = 0
        _tdly[:,:31] = 0; _tdly[:,-30:] = 0
        _tdly[:,77:93] = 0
        _ant = n.fft.ifft(_dly, axis=1)
        _tant = n.fft.ifft(_tdly, axis=1)
        _ant *= _wgt
        _tant *= _wgt
        #p.subplot(4, len(keys)+1, cnt+1)
        #waterfall(_ant, mode='lin')

        #p.subplot(4, len(keys)+1, 2*(len(keys)+1)+cnt+1)
        #waterfall(_tant, mode='lin')

        #p.subplot(4, len(keys)+1, len(keys)+1)
        _ant = _ant.sum(axis=0)
        _ant2 = n.abs(_ant)**2
        _wgt = _wgt.sum(axis=0)
        _wgt2 = n.abs(_wgt)**2
        step = []
        for thresh in range(1,_ant.size-1):
            left = n.sqrt(_ant2[:thresh].sum() / _wgt2[:thresh].sum().clip(1,n.Inf))
            right = n.sqrt(_ant2[thresh:].sum() / _wgt2[thresh:].sum().clip(1,n.Inf))
            step.append(left - right)
            if cnt < 4: step[-1] += .0005
        _avg = _ant / _wgt.clip(1,n.Inf)
        if _min is None: _min = _avg.copy()
        else: _min = n.where(n.abs(_avg) < n.abs(_min), _avg, _min)
        #p.plot(n.abs(_avg))#, label=filename)
        #p.plot(step)#, label='step'+filename)

        #p.subplot(4, len(keys)+1, 3*(len(keys)+1))
        _tant = _tant.sum(axis=0)
        _tant2 = n.abs(_tant)**2
        tstep = []
        for thresh in range(1,_tant.size-1):
            left = n.sqrt(_tant2[:thresh].sum() / _wgt2[:thresh].sum().clip(1,n.Inf))
            right = n.sqrt(_tant2[thresh:].sum() / _wgt2[thresh:].sum().clip(1,n.Inf))
            tstep.append(left - right - .0005)
        _tavg = _tant / _wgt.clip(1,n.Inf)
        #p.plot(n.abs(_tavg))#, label=filename)
        #p.plot(tstep)#, label='step'+filename)

        #p.subplot(4, len(keys)+1, 1*(len(keys)+1)+cnt+1)
        #waterfall(_dly, recenter=[1], drng=2)

        #p.subplot(4, len(keys)+1, 3*(len(keys)+1)+cnt+1)
        #waterfall(_tdly, recenter=[1], drng=2)
        ##waterfall(_dly[:,50:128], mode='phs')

        #p.subplot(4, len(keys)+1, 2*(len(keys)+1))
        ##p.plot(n.abs(n.fft.fft(_avg)[:128]), label=filename)
        #p.plot(n.abs(n.fft.fft(_avg)), label=filename)
        #p.ylim(0,.01)

        #p.subplot(4, len(keys)+1, 4*(len(keys)+1))
        #p.plot(n.abs(n.fft.fft(_tavg)), label=filename)

p.plot(_min)

#tot,totwgt = 0,0
#for bl in ant:
#    tot += ant[bl]
#    totwgt += wgt[bl]


