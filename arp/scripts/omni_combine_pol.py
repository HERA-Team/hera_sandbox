#! /usr/bin/env python
import numpy as n, sys

for fx in sys.argv[1:]:
    fy = fx.replace('xx','yy')
    print 'Reading', fx, fy
    nx, ny = n.load(fx), n.load(fy)
    d = {}
    d['lsts'] = nx['lsts']
    for k in nx:
        if not k.startswith('sep'): continue
        dx, dy = nx[k], ny[k]
        wx, wy = nx['wgt_'+k], ny['wgt_'+k]
        d[k] = .5 * (dx + dy)
        d['wgt_'+k] = .5 * (wx + wy)
    fI = fx.replace('xx','I')
    print 'Writing', fI
    n.savez(fI, **d)
