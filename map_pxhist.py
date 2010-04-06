#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

def rev(d):
    o = n.arange(d.size-1, -1, -1)
    return d[o]

for i,f in enumerate(sys.argv[1:]):
    m = a.map.Map(fromfits=f)

    wgts = n.where(m.wgt.map == 0, 1, m.wgt.map)
    flx = m.map.map / wgts
    wgt_sort = rev(n.argsort(m.wgt.map))
    wgt_sort = m.wgt.map[wgt_sort] / n.sum(m.wgt.map)
    csum_wgt = n.cumsum(wgt_sort)
    px = n.argwhere(csum_wgt > 1 - 1e-4)[0]
    thresh = wgt_sort[px]
    
    valid = n.where(m.wgt.map > thresh, 1, 0)
    flx = flx.compress(valid)

    #h, bins = n.histogram(n.log10(m.wgt.map.clip(1e-5,n.Inf)), bins=20)
    h, bins = n.histogram(n.log10(n.abs(flx)), bins=200)
    bins = .5 * (bins[:-1] + bins[1:])
    #order = n.argsort(flx)
    #flx = flx[order]
    p.subplot(len(sys.argv[1:]), 1, i+1)
    p.semilogy(bins, h, '.', label='src #')
    p.semilogy(bins, rev(n.cumsum(rev(h))), '-', label='srcs >')
    p.semilogy(bins, 10**bins * h, '.', label='flx #')
    p.semilogy(bins, n.cumsum(10**bins * h), '-', label='flx >')
    p.semilogy(bins, rev(n.cumsum(rev(10**bins * h))), '-', label='flx <')
    p.xlim(-4.5,4.5)
    p.ylim(1e-1,1e7)
    p.grid()
    p.legend(loc='lower left')
    p.ylabel('Counts')
p.xlabel('Log10(Jy)')
p.show()

