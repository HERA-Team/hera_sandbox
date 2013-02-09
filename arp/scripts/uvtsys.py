#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

freqs = aa.get_afreqs()
jy2T = C.pspec.jy2T(freqs)

_sum, _wgt = {}, {}
for filename in args:
    print 'Reading', filename,
    uv = a.miriad.UV(filename)
    print uv['sdf'], uv['inttime']
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    times = []
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            aa.set_jultime(t)
            times.append(t)
        bl = a.miriad.ij2bl(i,j)
        d[1:] = (d[:-1] - d[1:]) / n.sqrt(2)
        f[1:] = n.logical_or(f[:-1], f[1:])
        passband = aa.passband(i,j)
        d /= n.where(passband == 0, 1., passband)
        val = n.logical_not(f)
        #d *= jy2T * val * (2 * uv['sdf'] * 1e9 * uv['inttime'])**.5
        d *= jy2T * val
        _sum[bl] = _sum.get(bl,0) + n.abs(d)
        _wgt[bl] = _wgt.get(bl,0) + val

for bl in _sum:
    spec = _sum[bl] / _wgt[bl].clip(n.median(_wgt[bl])/2,n.Inf)
    #print a.miriad.bl2ij(bl), n.median(spec)
    p.semilogy(spec)
p.show()
