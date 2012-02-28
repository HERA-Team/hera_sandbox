#! /usr/bin/env python
"""
Filter visibilites per-baseline using a delay transform.
"""

import aipy as a, numpy as n, os, sys, optparse 

def gen_skypass_delay(aa, sdf, nchan, dw=0., max_bl_frac=1.):
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = max(dw, max_bl_frac * n.sqrt(n.dot(max_bl, max_bl)))
        dly,off = aa.get_phs_offset(i,j)[-2:]
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(n.round(uthresh)), int(n.round(lthresh))
        filters[bl] = (uthresh,lthresh)
        print (i,j), filters[bl]
    return filters

o = optparse.OptionParser()
o.set_usage('filter_sky.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=0,
    help='The number of delay bins to null. Default 0.')
o.add_option('--horizon', dest='horizon', type=float, default=0.,
    help='The additional scalar applied to the baseline length to determine horizon cutoff.  Default is 0 (no horizon filtering).')
o.add_option('--invert', dest='invert', action='store_true',
    help='Invert filter to preserve high delays rather than low.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
filters = gen_skypass_delay(aa, uv['sdf'], uv['nchan'], dw=opts.dw, max_bl_frac=opts.horizon)

for uvfile in args:
    uvofile = uvfile + 'F'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    def mfunc(uv, p, d, f):
        crd,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        val = n.logical_not(f).astype(n.float)
        gain = n.sqrt(n.average(val**2))
        if n.all(val == 0): return p, d, f
        ker = n.fft.ifft(val)
        _d = n.fft.ifft(d)
        _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        _d += info['res'] / gain
        uthresh,lthresh = filters[bl]
        if opts.invert:
            if lthresh < 0:
                _d[:uthresh] = 0
                _d[lthresh:] = 0
            else:
                _d[lthresh:uthresh] = 0
        else:
            if lthresh < 0:
                _d[uthresh:lthresh] = 0
            else:
                _d[uthresh:] = 0
                _d[:lthresh] = 0
        d = n.fft.fft(_d) * val
        return p, d, f

    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    # Apply the pipe to the data
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')

