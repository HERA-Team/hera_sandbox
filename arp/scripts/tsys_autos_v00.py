#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, max=True, drng=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

data = {}
curtime = None
for filename in args:
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d in uv.all():
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
        #d /= aa.passband(i,j)
        d = n.abs(d)
        if i == j:
            try: data[bl].append(d)
            except(KeyError): data[bl] = [d]
        else:
            try: data[bl].append(d)
            except(KeyError): data[bl] = [d]


print len(data)
for cnt,bl in enumerate(data):
    p.subplot(len(data), 1, cnt+1)
    i,j = a.miriad.bl2ij(bl)
    d = n.ma.array(data[bl])
    dd = d[1:-1] - .5 * (d[:-2] + d[2:])
    #rms = n.ma.std(dd, axis=0)
    if i != j:
        bli = a.miriad.ij2bl(i,i)
        blj = a.miriad.ij2bl(j,j)
        d = n.sqrt(n.ma.array(data[bli]) * n.ma.array(data[blj]))
    avg = n.ma.average(d[1:-1], axis=0)
    #for cnt in range(0):
    #    dd = n.ma.masked_where(n.abs(dd) > 2*rms, dd)
    #    d[1:-1] = n.ma.masked_where(n.abs(dd) > 2*rms, d[1:-1])
    #    rms = n.ma.std(dd, axis=0)
    #    avg = n.ma.average(d[1:-1], axis=0)
    if False:
        scale = n.sqrt(2 * uv['sdf'] * 1e9 * uv['inttime']) / n.sqrt(3./2)
        p.plot(rms / avg * scale)
    else:
        #d = n.log10(n.abs(d))
        d = n.log10(n.abs(dd) / d[1:-1] * n.sqrt(2 * uv['sdf'] * 1e9 * uv['inttime']) / n.sqrt(3./2))
        if opts.max is None: mx = d.max()
        else: mx = opts.max
        if opts.drng is None: mn = max(d.min(), -10)
        else: mn = mx - opts.drng
        p.imshow(d, vmax=mx, vmin=mn, aspect='auto')
        p.colorbar(shrink=.5)
p.show()
