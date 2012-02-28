#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, chan=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.select_chans(chans)
del(uv)

bl_gom = a.miriad.ij2bl(1,1)
dat = {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)), d, f in uv.all(raw=True):
        d = d.take(chans)
        f = f.take(chans)
        bl = a.miriad.ij2bl(i,j)
        if not dat.has_key(bl): dat[bl] = []
        if bl == bl_gom:
            dat[bl].append(d)
        else:
            dat[bl].append(n.where(f, 0, d))

gom = n.abs(n.array(dat[bl_gom]))
del(dat[bl_gom])
for bl in dat:
    i,j = a.miriad.bl2ij(bl)
    gain = aa.passband(i,j)
    dat[bl] = n.array(dat[bl]) / gain
    #p.semilogy(n.average(n.abs(dat[bl]) / gom, axis=1))
    p.semilogy(n.average(n.abs(dat[bl]), axis=1))
p.semilogy(1e4*n.average(gom, axis=1))
p.show()
