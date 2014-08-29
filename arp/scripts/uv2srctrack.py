#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True, cal=True, chan=True)
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
src = cat.values()[0]

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
fq = n.average(aa.get_afreqs().take(chans))
del(uv)

times = []
x,y,z = [],[],[]
dat,wgt = [], []

for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    curtime = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            aa.set_jultime(t)
            src.compute(aa)
            xi,yi,zi = src.get_crds('top', ncrd=3)
        w = n.logical_not(f.take(chans))
        d = d.take(chans)
        times.append(t)
        x.append(xi); y.append(yi); z.append(zi)
        dat.append(n.sum(d))
        wgt.append(n.sum(w))

times = n.array(times)
x,y,z = n.array(x), n.array(y), n.array(z)
dat, wgt = n.array(dat), n.array(wgt)

#filename = 'srctrack__%s.npz' % src.src_name
filename = '%s__srctrack.npz' % src.src_name
print 'Writing', filename
n.savez(filename, times=times, x=x, y=y, z=z, spec=dat, wgt=wgt, freq=fq)
