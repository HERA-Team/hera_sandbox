#! /usr/bin/env python
import numpy as n, pylab as p, capo as C
import sys

for f in sys.argv[1:]:
    meta,gains,vismdl,xtalk = C.omni.from_npz(f)
    p.subplot(131)
    nchan = meta['chisq'].shape[1]
    C.arp.waterfall(meta['chisq'], drng=3); p.colorbar()

    nant = max(gains.keys()) + 1
    amps = n.zeros((nant,nchan), dtype=n.complex)
    vars = n.zeros((nant,nchan), dtype=n.complex)
    for i in gains:
        amps[i] = n.median(gains[i], axis=0)
        vars[i] = n.var(gains[i], axis=0)
    p.subplot(132)
    p.title(f)
    C.arp.waterfall(amps, drng=2); p.colorbar()
    p.subplot(133)
    C.arp.waterfall(vars, drng=4); p.colorbar()
    p.show()

