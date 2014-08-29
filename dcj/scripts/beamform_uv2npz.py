#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, re

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, dec=True, src=True, chan=True)
o.add_option('--addsrc', dest='addsrc', action='store_true',
    help='For beamform files that are residuals, add back in the source model.')
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
cat.compute(aa)
src = cat.values()[0]
print src

del(uv)

# Need to stack: spec, times, wgts, x, y, z, path name
# File name is srcname__*

times = []
spec,wgts = [],[]
x,y,z = [],[],[]
beam = []
lst = []
flux = []
flags = []
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
            src.compute(aa)
            xi,yi,zi = src.get_crds('top')
            print t,aa.sidereal_time()
        if src.alt < opts.altmin * a.img.deg2rad: continue
        x.append(xi); y.append(yi); z.append(zi)
        times.append(t)
        d,w = d.take(chans), n.logical_not(f.take(chans)).astype(n.int)
        bm = aa[0].bm_response((xi,yi,zi), pol=opts.pol[0]) * aa[0].bm_response((xi,yi,zi), pol=opts.pol[1])
        lst.append(uv['lst'])
        if opts.addsrc:     
            d += src.jys.flatten() * bm.flatten()
        flux.append(n.sum(d))
        spec.append(d.flatten())
        wgts.append(n.sum(w))
        beam.append(bm.flatten())
        flags.append(f.take(chans).astype(n.int))

x,y,z = n.array(x), n.array(y), n.array(z)
spec,wgts = n.array(spec), n.array(wgts)
times = n.array(times)
beam =n.array(beam)
flags = n.array(flags)
afreqs = aa.get_afreqs()
#npzfile = src.src_name + '_spec.npz'
npzfile = src.src_name+'__srctrack_%.0f-%.0f.npz'%(n.min(times),n.max(times))
print 'Writing spectrum to', npzfile
n.savez(npzfile, spec=spec, flux=flux, freq=afreqs, times=times, x=x, y=y, z=z, wgts=wgts,beam=beam,lst=lst,flags=flags)
