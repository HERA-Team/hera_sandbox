#! /usr/bin/env python
import numpy as n, pylab as p, capo as C
import sys

pol = 'xx'
for f in sys.argv[1:]:
    meta,gains,vismdl,xtalk,jds,lsts,freqs = C.omni.from_npz(f)
    p.subplot(131)
    nchan = meta[pol]['chisq'].shape[1]
    C.arp.waterfall(meta[pol]['chisq'], drng=3); p.colorbar()

    nant = max(gains[pol[0]].keys()) + 1
    amps = n.zeros((nant,nchan), dtype=n.complex)
    chis = n.zeros((nant,nchan), dtype=n.complex)
    for i in gains[pol[0]]:
        amps[i] = n.abs(n.median(gains[pol[0]][i], axis=0))
        chis[i] = n.median(meta[pol]['chisq%d' % i], axis=0)
    p.subplot(132)
    p.title(f)
    C.arp.waterfall(amps, mode='lin', mx=1.5, drng=1); p.colorbar()
    p.subplot(133)
    C.arp.waterfall(chis / n.where(amps > 0, amps, 1), drng=4); p.colorbar()
    p.show()

