#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, dec=True, 
    src=True, chan=True)
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
if opts.chan is None: opts.chan = 'all'
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
srclist,cutoff,catalogs, = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
src = cat.values()[0]
del(uv)

spec, swgt = 0, 0
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    curtime = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            cat.compute(aa)
            if src.alt < opts.altmin * a.img.deg2rad: continue
            s_eq = cat.get_crds('eq', ncrd=3)
            aa.sim_cache(s_eq)
        if src.alt < opts.altmin * a.img.deg2rad: continue

        d,f = d.take(chans), f.take(chans)
        wgt = aa.bm_response(i,j,pol=opts.pol).squeeze()
        d = n.where(f, 0, d); wgt = n.where(f, 0, wgt)
        # Optimal SNR: down-weight beam-attenuated data 
        # by another factor of the beam response.
        d *= wgt; wgt *= wgt
        spec += d; swgt += wgt

#spec = spec.real / swgt
spec = n.abs(spec) / swgt
valid = n.logical_not(n.isnan(spec))
spec = spec.compress(valid)
afreqs = aa.get_afreqs().compress(valid)
npzfile = src.src_name + '_spec.npz'
print 'Writing spectrum to', npzfile
n.savez(npzfile, spec=spec, freq=afreqs)
