#! /usr/bin/env python
import numpy as n, pylab as p, capo as C
import sys

for f in sys.argv[1:]:
    npz = n.load(f)
    p.subplot(131)
    nchan = npz['xx,chisq'].shape[1]
    C.arp.waterfall(npz['xx,chisq'], drng=3); p.colorbar()

    gains = [k for k in npz.files if k.find('gains') != -1]
    ants = [int(k.split(',')[-1]) for k in gains]
    nant = max(ants) + 1
    amps = n.zeros((nant,nchan), dtype=n.complex)
    vars = n.zeros((nant,nchan), dtype=n.complex)
    for g in gains:
        pol,_,i = g.split(',')
        lbl,i = i+pol, int(i)
        amps[i] = n.median(npz[g], axis=0)
        vars[i] = n.var(npz[g], axis=0)
    p.subplot(132)
    p.title(f)
    C.arp.waterfall(amps, drng=2); p.colorbar()
    p.subplot(133)
    C.arp.waterfall(vars, drng=4); p.colorbar()
    p.show()

