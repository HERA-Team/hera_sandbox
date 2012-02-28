#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, chan=True)
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
sdf = uv['sdf']
inttime = uv['inttime']
del(uv)


dat,val = {}, {}
for uvfile in args:
    print 'Reading', uvfile
    uvi = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uvi, opts.ant)
    # Gather all data and each time step
    for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        pol = uvi['pol']
        if not dat.has_key(pol): dat[pol], val[pol] = {}, {}
        if not dat[pol].has_key(bl): dat[pol][bl], val[pol][bl] = [], []
        val[pol][bl].append(n.logical_not(f.take(chans)).astype(n.int))
        dat[pol][bl].append(d.take(chans))

# Generate statistical mask
for pol in dat:
  for bl in dat[pol]:
    d,v = n.array(dat[pol][bl]), n.array(val[pol][bl])
    sum2_f = (n.abs(d)**2).sum(axis=0)
    wgt2_f = (n.abs(v)**2).sum(axis=0)
    rms_f = n.sqrt(sum2_f / wgt2_f.clip(1,n.Inf))
    for i,color in zip([256,64,16,4,1],'rmbck'):
        dim = d.shape[0] / i
        sum_f, wgt_f = 0, 0
        for j in range(0, d.shape[0], dim):
            dj,vj = d[j:j+dim], v[j:j+dim]
            sum_f += n.abs(dj.sum(axis=0))**2
            wgt_f += n.abs(vj.sum(axis=0))**2
        avg_f = n.sqrt(sum_f / n.where(wgt_f > 0, wgt_f, 1))
        p.plot(rms_f / n.sqrt(n.float(dim)), color+':')
        p.semilogy(n.abs(avg_f), color+'-')
    p.show()
    d = n.array(dat[pol][bl]).flatten()
    v = n.array(val[pol][bl]).flatten()
    i, j = a.miriad.bl2ij(bl)
    d = n.abs(d.compress(v))
    hist,bins = n.histogram(d * n.sqrt(sdf * 1e9 * inttime) * 2e-3, bins=50)
    p.plot(bins[:-1], hist)
    p.xlabel('Temperature [K]')
    p.ylabel('Counts')
    p.show()



