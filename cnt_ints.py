#! /usr/bin/env python
"""Count the number of unflagged integrations in each channel range."""
import aipy as a, numpy as n, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, chan=True, cal=True, src=True, 
    dec=True, pol=True)
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
o.add_option('--no_bm', dest='no_bm', action='store_true',
    help="Don't weight integrations by beam response.")
o.add_option('--samples', dest='samples', action='store_true',
    help="Count samples (i.e. each baseline seperately) instead of integrations.")
opts,args = o.parse_args(sys.argv[1:])
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'], concat=False)
del(uv)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
src = cat.values()[0]
ints = 0
curtime = None
for filename in args:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            cat.compute(aa)
            if src.alt < opts.altmin * a.img.deg2rad: continue
            s_eq = cat.get_crds('eq', ncrd=3)
            aa.sim_cache(s_eq)
            if opts.no_bm: wgt = 1
            else: wgt = aa.bm_response(0,0, pol=opts.pol).squeeze()
            if not opts.samples: ints += n.logical_not(f).astype(n.int) * wgt**2
        if opts.samples: ints += n.logical_not(f).astype(n.int) * wgt**2
print
print "'%s' : {" % opts.src
for ch in chans:
    print '%d : %10.1f,' % (ch[0], ints.take(ch).sum() / float(len(ch)))
print '},'
