#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, dec=True, 
    src=True, chan=True)
o.add_option('-d','--deg', dest='deg', type='int', default=8,
    help='Degree of polynomial to fit to bandpass function.')
o.add_option('-q','--quiet', dest='quiet', action='store_true',
    help='Do not plot anything (be visually quiet).')
o.add_option('-o','--outfile', dest='outfile', 
    help='Output file to save plot to.  If none, opens window.')
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
o.add_option('--minuv', dest='minuv', type='float', default=0,
    help='Minimum uv length (in wavelengths) for a baseline to be included.')
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
        u,v,w = aa.gen_uvw(i, j, src)
        tooshort = n.where(n.sqrt(u**2+v**2) < opts.minuv, 1, 0).squeeze()
        if n.all(tooshort): continue

        d,f = d.take(chans), f.take(chans)
        f = n.logical_or(f, tooshort)
        d /= aa.passband(i,j)
        # If source is resolved, downweight baselines that resolve it
        try:
            res = n.abs(aa.gen_phs(src, i, j, ionref=src.ionref, 
                srcshape=src.srcshape, resolve_src=True))
        except(a.phs.PointingError): continue
        d = aa.phs2src(d, src, i, j)
        wgt = res * aa.bm_response(i,j,pol=opts.pol).squeeze()
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
#vspec = spec
#vafreqs = afreqs
dspec = spec - a.rfi.remove_spikes(spec, order=8)
sig = n.std(dspec)
valid = n.where(n.abs(dspec) < 2*sig, 1, 0)
vspec = spec.compress(valid)
vafreqs = afreqs.compress(valid)
src.update_jys(vafreqs)
src_poly = n.polyfit(n.log10(vafreqs/src.mfreq), n.log10(vspec), deg=1)
bp = n.sqrt(vspec / src.jys)
bp_poly = n.polyfit(vafreqs, bp, deg=opts.deg)
bp_fit = n.polyval(bp_poly, vafreqs).clip(.5,2)**2
n.set_printoptions(threshold=n.nan)
#print args
print 'bp =', list(bp_poly)
print "'%s':" % opts.src + "{ 'jys':10**%f, 'index':  %f , }," % (src_poly[-1], src_poly[-2])
print 'RMS residual:', n.sqrt(n.average((vspec - 10**n.polyval(src_poly, n.log10(vafreqs/src.mfreq)))**2))

if not opts.quiet:
    if not opts.outfile is None:
        import matplotlib
        matplotlib.use('Agg')
    import pylab as p
    p.loglog(afreqs, spec, 'k^', label='Measured')
    p.loglog(vafreqs, src.jys, 'b:', 
        label='%f, %s' % (src._jys, str(src.index)))
    p.loglog(vafreqs, vspec / bp_fit, 'g.', label='Post-Fit')
    p.loglog(afreqs, 10**n.polyval(src_poly, n.log10(afreqs/src.mfreq)), 'r:', 
        label='Fit Power Law')
    p.xticks(n.arange(.1,.2,.01))
    p.xlim(afreqs[0], afreqs[-1])
    p.ylim(1,2e4)
    p.grid()
    p.title(src.src_name)
    #p.legend()
    if not opts.outfile is None: p.savefig(opts.outfile)
    else: p.show()
