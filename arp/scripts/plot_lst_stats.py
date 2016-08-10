#! /usr/bin/env python

import aipy as a, capo as C, pylab as p, numpy as n
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

cnt, var = {}, {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        bl = '%d,%d,%d' % (i,j,uv['pol'])
        cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
        var[bl] = var.get(bl, []) + [uv['var']]

for bl in cnt:
    cnt[bl] = n.array(cnt[bl])
    var[bl] = n.array(var[bl])

for bl in cnt:
    print bl
    p.subplot(121); C.arp.waterfall(cnt[bl], mode='lin'); p.colorbar(shrink=.5)
    p.subplot(122); C.arp.waterfall(var[bl], mode='log'); p.colorbar(shrink=.5)
    p.show()
