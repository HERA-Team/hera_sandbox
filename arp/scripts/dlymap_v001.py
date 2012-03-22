#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, chan=True, dec=True)
o.add_option('--mfreq', dest='mfreq', type='float', default=.150,
    help='Frequency to use for sub-bin delay extrapolations.  Default .150 GHz')

opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
(crd,t,(i,j)),d,f = uv.read(raw=True)
#chans = a.scripting.parse_chans(opts.chan uv['nchan'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.set_jultime(t)
#aa.select_chans(chans)
afreqs = aa.get_afreqs()
window = a.dsp.gen_window(afreqs.size, window='blackman-harris')
bins = n.fft.fftfreq(afreqs.size, afreqs[1] - afreqs[0])
dbin = bins[1] - bins[0]

h = a.map.Map(nside=512)
px = n.arange(h.npix())
eq_J2000 = n.array(h.px2crd(px))
m_precess = a.coord.convert_m('eq','eq',oepoch=aa.epoch)
eq = n.dot(m_precess, eq_J2000)

import pylab as p

for filename in args:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    curtime = None
    for (crd,t,(i,j)), d, f in uv.all(raw=True):
        print i,j, t, curtime
        if t != curtime:
            aa.set_jultime(t)
            m = n.dot(aa.eq2top_m, m_precess)
            tx,ty,tz = n.dot(aa.eq2top_m,eq)
            valid = n.where(tz > 0)
            tx,ty,tz = tx[valid], ty[valid], tz[valid]
            curtime = t
        d = n.where(f, 0, aa.phs2src(d, 'z', i, j))
        d = n.conj(d)
        #d = n.where(f, 0, n.ones_like(d))
        if n.all(f): continue
        w = n.logical_not(f).astype(n.float)
        _d, _w = n.fft.ifft(d*window), n.fft.ifft(w*window)
        if True:
            _d, info = a.deconv.clean(_d, _w, tol=1e-5)
            gain = n.sqrt(n.average(w))
            _d += info['res'] / gain
        _d = n.real(_d)
        bx,by,bz = aa.get_baseline(i,j,src='e')
        tau = bx*eq[0][valid] + by*eq[1][valid] + bz*eq[2][valid]
        taubin = n.around(tau / dbin).astype(n.int)
        mdat = n.abs(_d[taubin])
        #mdat = tau
        h.add(px[valid], n.ones_like(mdat), mdat)
        
h.to_fits('out.fits', clobber=True)
