#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, src=True, dec=True)
o.add_option('-b', '--beam', action='store_true',
    help='Normalize by the primary beam response in the direction of the specified source.')
o.add_option('-f', '--srcflux', action='store_true',
    help='Normalize by the spectrum of the specified source.')
o.add_option('--minuv', type='float', default=0,
    help='Minimum uv length (in wavelengths) for a baseline to be included.')
o.add_option('--maxuv', type='float', default=n.Inf,
    help='Maximum uv length (in wavelengths) for a baseline to be included.')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
src = cat.values()[0]
del(uv)

for filename in args:
    print filename, '->', filename+'.bm_'+src.src_name
    if os.path.exists(filename+'.bm_'+src.src_name):
        print '    File exists, skipping.'
        continue
    dbuf,cbuf = {}, {}
    curtime = None
    print '    Summing baselines...'
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants, opts.pol)
    aa.set_active_pol(opts.pol)
    uvi.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
            src.compute(aa)
            if opts.srcflux: src_spec = src.get_jys()
            if opts.beam: 
                s_eq = cat.get_crds('eq', ncrd=3)
                aa.sim_cache(s_eq)
        try:
            d = aa.phs2src(d, src, i, j)
            u,v,w = aa.gen_uvw(i, j, src)
            tooshort = n.where(n.sqrt(u**2+v**2) < opts.minuv, 1, 0).squeeze()
            toolong = n.where(n.sqrt(u**2+v**2) > opts.maxuv, 1, 0).squeeze()
            dont_use = n.logical_or(tooshort, toolong)
            if n.all(dont_use):
                #print i,j, 'too short:',
                #print n.average(n.sqrt(u**2+v**2)), '<', opts.minuv
                continue
            f = n.logical_or(f, dont_use)
            gain = aa.passband(i,j)
            if opts.beam: gain *= aa.bm_response(i,j,pol=opts.pol).squeeze()
            if opts.srcflux: gain *= src_spec
            d /= gain
        except(a.phs.PointingError): d *= 0
        dbuf[t] = dbuf.get(t, 0) + n.where(f, 0, d)
        cbuf[t] = cbuf.get(t, 0) + n.logical_not(f).astype(n.int)
    uvi.rewind()

    print '    Writing output file'
    curtime = None
    def mfunc(uv, p, d, f):
        global curtime
        uvw,t,(i,j) = p
        if t != curtime:
            curtime = t
            cnt = cbuf[t].clip(1,n.Inf)
            f = n.where(cbuf[t] == 0, 1, 0)
            d = dbuf[t] / cnt
            return (uvw,t,(i,j)), d, f
        else: return (uvw,t,(1,1)), None, None
        
    uvo = a.miriad.UV(filename+'.bm_'+src.src_name, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='BEAMFORM: src=%s ant=%s srcflux=%s minuv=%s\n' % (opts.src, opts.ant, opts.srcflux, opts.minuv))
    del(uvi); del(uvo)
    
    
