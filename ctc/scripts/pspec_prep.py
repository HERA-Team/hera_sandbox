#! /usr/bin/env python
"""
Filter visibilites per-baseline using a delay transform.
"""

import aipy as a, numpy as n, os, sys, optparse 
import capo as C

def gen_skypass_delay(aa, sdf, nchan, max_bl_add=0.):
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        #print (i,j), max_bl,
        max_bl = max_bl_add + n.sqrt(n.dot(max_bl, max_bl))
        uthresh, lthresh = max_bl/bin_dly + 1.5, -max_bl/bin_dly - 0.5
        uthresh, lthresh = int(n.ceil(uthresh)), int(n.floor(lthresh))
        #print max_bl, uthresh, lthresh
        filters[bl] = (uthresh,lthresh)
        #if i == 0 and j == 60: print (i,j), max_bl / max_bl_frac, max_bl_frac, filters[bl]
    return filters

o = optparse.OptionParser()
o.set_usage('pspec_prep.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--nophs', dest='nophs', action='store_true',
    help='Do not phase to zenith bin.')
o.add_option('--nogain', dest='nogain', action='store_true',
    help='Do not normalize gain.')
o.add_option('--window', dest='window', default='none',
    help='DSP window to use.  Default: none')
o.add_option('--windowchan', dest='windowchan', default=None,
    help='Frequency channel range (inclusive) for the window function (ex: 16_182). Default uses the entire frequency range.')
o.add_option('--horizon', dest='horizon', type=float, default=0.,
    help='An additional additive term (in ns) applied to the baseline length to determine horizon cutoff.  Default is 0.')
o.add_option('--clean', dest='clean', type='float', default=1e-5,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination.  Default 1e-5')
o.add_option('--model', dest='model', action='store_true',
    help='Return the foreground model summed with the residuals (in Fourier space).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
filters = gen_skypass_delay(aa, uv['sdf'], uv['nchan'], max_bl_add=opts.horizon)

for uvfile in args:
    if opts.model: uvofile = uvfile + 'F'
    else: uvofile = uvfile + 'B'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.add_var('bin','d')

    if opts.windowchan != None: #smaller window function range
        f1,f2 = opts.windowchan.split('_')[0],opts.windowchan.split('_')[1]
        frange = int(f2)-int(f1)+1
        window_short = a.dsp.gen_window(frange, window=opts.window)
        end1 = n.zeros(int(f1))
        end2 = n.zeros(uvi['nchan']-int(f2)-1)
        window = n.concatenate([end1,window_short,end2])
    else:
        window = a.dsp.gen_window(uvi['nchan'], window=opts.window)
    
    curtime, zen = None, None
    def mfunc(uv, p, d, f):
        global curtime,zen
        crd,t,(i,j) = p
        if t != curtime:
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            ubin,vbin,lstbin = C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst))
            zen = a.phs.RadioFixedBody(lstbin, aa.lat)
            zen.compute(aa)
            curtime = t
            print t
        pol = a.miriad.pol2str[uv['pol']]
        aa.set_active_pol(pol)
        u,v,w = aa.gen_uvw(i,j, src=zen)
        u,v = u.flatten()[-1], v.flatten()[-1]
        conj = False
        if u < 0: 
            u,v = -u, -v
            if not opts.nophs: conj = True
        uvo['bin'] = n.float(C.pspec.uv2bin(u, v, aa.sidereal_time()))
        if i == j: return p, d, f
        bl = a.miriad.ij2bl(i,j)
        w = n.logical_not(f).astype(n.float)
        if n.average(w) < .5: return p, n.zeros_like(d), n.ones_like(f)
        if not opts.nophs: d = aa.phs2src(d, zen, i, j)
        if not opts.nogain: d /= aa.passband(i,j)
        if conj: d = n.conj(d)
        d *= w
        _d = n.fft.ifft(d * window)
        _w = n.fft.ifft(w * window)
        uthresh,lthresh = filters[bl]
        area = n.ones(_d.size, dtype=n.int)
        area[uthresh:lthresh] = 0
        _d_cl, info = a.deconv.clean(_d, _w, tol=opts.clean, area=area, stop_if_div=False, maxiter=100)
        if opts.model:
            d_mdl = n.fft.fft(_d_cl + info['res'])
            f = n.zeros_like(d_mdl)
            return p, d_mdl, f
        else:
            d_mdl = n.fft.fft(_d_cl)
            d_res = d - d_mdl * w
            return p, d_res, f
 
    # Apply the pipe to the data
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
