#! /usr/bin/env python
"""Count the number of unflagged integrations in each channel range."""
import aipy as a, numpy as n, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True, cal=True, src=True, dec=True)
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
opts,args = o.parse_args(sys.argv[1:])
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'], concat=False)
del(uv)
srclist,cutoff = a.scripting.parse_srcs(opts.src)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff)
src = cat.values()[0]
ints = 0
curtime = None
for filename in args:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = a.miriad.UV(filename)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            src.compute(aa)
            if src.alt < opts.altmin * a.img.deg2rad: continue
            ints += n.logical_not(f).astype(n.int)
print
print "'%s' : {" % opts.src
for ch in chans:
    print '%d : %10.1f,' % (ch[0], ints.take(ch).sum() / float(len(ch)))
print '},'
