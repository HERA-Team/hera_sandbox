#! /usr/bin/env python
import aipy as a, numpy as n
#from arp_paper import *

aspec = {}
info = {}
for filename in ['zen.2455332.81387.uvcAR','zen.2455332.54941.uvcAR', 'zen.2455332.31975.uvcAR']:
    aspec[filename], info[filename] = {}, {}
    uv = a.miriad.UV(filename)
    chan = n.arange(uv['nchan'])
    freq = chan * uv['sdf'] + uv['sfreq']
    times = []
    a.scripting.uv_selector(uv, '1_1', 'xx')
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]: times.append(t)
        bl = a.miriad.ij2bl(i,j)
        v = n.logical_not(f).astype(n.float)
        _sum,_wgt = aspec[filename].get(i,(0,0))
        _sum += d * v
        _wgt += v
        aspec[filename][i] = (_sum, _wgt)
    info[filename]['freq'] = freq

import pylab as p
for filename in aspec:
    for i in aspec[filename]:
        _sum,_wgt = aspec[filename][i]
        _avg = n.abs(_sum / _wgt.clip(1,n.Inf))
        _avg /= n.average(_avg[140:150])
        aspec[filename][i] = _avg

for filename in aspec:
    for i in aspec[filename]:
        #p.semilogy(info[filename]['freq'], n.abs(aspec[filename][i] - aspec['zen.2455332.81387.uvcAR'][i]) / aspec['zen.2455332.81387.uvcAR'][i])
        #p.semilogy(info[filename]['freq'], n.abs(aspec[filename][i]), label=filename+'_%d'%i)
        p.semilogy(n.abs(aspec[filename][i]), label=filename+'_%d'%i)

p.legend()
p.grid()
p.show()
