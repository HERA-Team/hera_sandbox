#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, os
#import pylab as p

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print 'Reading', filename
    buf = {}
    uv = a.miriad.UV(filename)
    buf['freqs'] = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    a.scripting.uv_selector(uv, ants=opts.ant)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        bl = str(a.miriad.ij2bl(i,j))
        dbl, wbl = 'd_'+bl, 'w_'+bl
        if not buf.has_key(dbl): buf[dbl],buf[wbl] = [],[]
        buf[dbl].append(n.abs(d))
        buf[wbl].append(n.logical_not(f))
    for dbl in buf:
        if not dbl.startswith('d_'): continue
        wbl = 'w'+dbl[1:]
        buf[wbl] = n.array(buf[wbl])
        buf[dbl] = n.array(buf[dbl]) * buf[wbl]
        wgt = buf[wbl].sum(axis=0).clip(1,n.Inf)
        avg = buf[dbl].sum(axis=0) / wgt
        avg_ = n.reshape(avg, (1,avg.size))
        dif = (buf[dbl] - avg_) * buf[wbl]
        std = n.sqrt(n.sum(dif**2, axis=0) / wgt)
        valid = n.where(n.abs(dif) > 2*std, 0, 1)
        buf[dbl] *= valid
        buf[wbl] *= valid
        buf[wbl] = buf[wbl].sum(axis=0)
        buf[dbl] = buf[dbl].sum(axis=0)
        #p.plot(buf[dbl] / buf[wbl].clip(1,n.Inf))
    print '    ->', filename+'__autoex.npz'
    n.savez(filename+'__autoex.npz', **buf)
#p.show()
