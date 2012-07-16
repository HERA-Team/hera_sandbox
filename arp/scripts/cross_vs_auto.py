#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys

times,d,f = C.arp.get_dict_of_uv_data(sys.argv[1:], 'all', 'xx', verbose=True)

for bl in d:
    i,j = a.miriad.bl2ij(bl)
    if i == j: continue
    dat = d[bl]['xx'] / n.sqrt(d[a.miriad.ij2bl(i,i)]['xx'] * d[a.miriad.ij2bl(j,j)]['xx'])
    C.arp.waterfall(dat, mx=-1.5, drng=2)
    p.colorbar(shrink=.5)
    p.show()
