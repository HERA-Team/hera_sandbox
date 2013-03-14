#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, os
import scipy.interpolate

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, src=True, dec=True)
o.add_option('--thumb', 
    help='Thumbnail to use as the source model')
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

npz = n.load(opts.thumb)
srcmdl = npz['img']
dra,ddec = npz['dra_ddec']
srcuv = n.fft.fft2(srcmdl); srcuv /= n.abs(srcuv).max()
u = n.fft.fftfreq(srcmdl.shape[0], dra)
v = n.fft.fftfreq(srcmdl.shape[1], ddec)
srcuv = n.abs(srcuv)
print srcuv.shape, u.shape, v.shape
print 'Creating interpolation'
uv_resp = scipy.interpolate.RectBivariateSpline(n.fft.fftshift(u),n.fft.fftshift(v),n.fft.fftshift(srcuv))
print 'Done creating interpolation'


for filename in args:
    print filename, '->', filename+'.bm_'+src.src_name
    if os.path.exists(filename+'.bm_'+src.src_name):
        print '    File exists, skipping.'
        continue
    dbuf,wbuf = {}, {}
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
        try:
            d = aa.phs2src(d, src, i, j)
            u,v,w = aa.gen_uvw(i, j, src)
            tooshort = n.where(n.sqrt(u**2+v**2) < opts.minuv, 1, 0).squeeze()
            toolong = n.where(n.sqrt(u**2+v**2) > opts.maxuv, 1, 0).squeeze()
            dont_use = n.logical_or(tooshort, toolong)
            if n.all(dont_use): continue
            u,v = u.flatten(), v.flatten()
            w = n.array([uv_resp(u0,v0) for u0,v0 in zip(u,v)]).flatten()
            #print u[u.size/2],v[u.size/2],w[u.size/2]
            w = n.where(dont_use, 0, w*n.logical_not(f))
            gain = aa.passband(i,j)
            d /= gain
        except(a.phs.PointingError): w = n.zeros_like(f)
        dbuf[t] = dbuf.get(t, 0) + d * n.conj(w)
        wbuf[t] = wbuf.get(t, 0) + w * n.conj(w)
    uvi.rewind()

    print '    Writing output file'
    curtime = None
    def mfunc(uv, p, d, f):
        global curtime
        uvw,t,(i,j) = p
        if t != curtime:
            curtime = t
            wgt = wbuf[t]
            f = n.where(wbuf[t] == 0, 1, 0)
            d = n.where(f, 0, dbuf[t] / wgt)
            return (uvw,t,(i,j)), d, f
        else: return (uvw,t,(1,1)), None, None
        
    uvo = a.miriad.UV(filename+'.bm_'+src.src_name, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='BEAMFORM: ' + ' '.join(sys.argv))
    del(uvi); del(uvo)
    
    
