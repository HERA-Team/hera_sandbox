#! /usr/bin/env python
"""
Filter in fringe-rate to select (or de-select) fringe rates that correspond to sources fixed
to the celestial sphere.
Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse
import capo as C
import pylab

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--maxbl', type='float', default=300.,
    help='Maximum baseline length (in meters) for purpose of calculating maximum fringe rates. Default 300.')
o.add_option('--lat', type='float', default=-30.,
    help='Latitude of array in degrees.  Default -30.')
o.add_option('--rmsky', action='store_true',
    help='Instead of retaining the data corresponding to the sky, remove it.')
o.set_usage('fringe_rate_filter_coarse.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
inttime = uv['inttime']
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

def sky_fng_thresh(bl_len, inttime, nints, freq, lat):
    '''For bl_ew_len (the east/west projection) in ns, return the (upper,negative,lower) fringe rate bins 
    that geometrically correspond to the sky.'''
    bin_fr = 1. / (inttime * nints)
    max_fr = freq * bl_len * 2*n.pi / a.const.sidereal_day
    neg_fr = max_fr * n.cos(n.pi/2 + n.abs(lat))
    nthr = neg_fr / bin_fr
    uthr = max_fr / bin_fr
    uthr, nthr = n.ceil(uthr).astype(n.int), int(n.floor(nthr))
    return uthr, nthr

def triplets(seq):
    for i in range(1,len(seq)-1): yield seq[i-1:i+2]
    return

data, wgts = {}, {}

for files in triplets(args):
    times, dat, flg = C.arp.get_dict_of_uv_data(files, opts.ant, opts.pol, verbose=True)
    # Variables: ufr (upper fringe rate), nfr (negative fringe rate)
    ufr,nfr = sky_fng_thresh(opts.maxbl, inttime, len(times), fqs.max(), opts.lat*a.img.deg2rad)

    window = a.dsp.gen_window(times.shape[0], 'blackman-harris')
    #window = 1
    area = n.ones(times.shape, dtype=n.int)
    area[ufr+1:nfr] = 0
    for bl in dat:
        if not data.has_key(bl): data[bl],wgts[bl] = {}, {}
        for pol in dat[bl]:
            if not data[bl].has_key(pol): data[bl][pol],wgts[bl][pol] = {}, {}
            d = n.where(flg[bl][pol], 0, dat[bl][pol])
            #w = n.logical_not(flg[bl][pol]).astype(n.float)
            w = n.where(n.abs(d) == 0, 0., 1.)
            for ch in range(d.shape[1]):
                _d,_w = n.fft.ifft(d[:,ch]*window), n.fft.ifft(w[:,ch]*window)
                gain = n.abs(_w[0])
                #print ch, nfr, ufr, gain
                if gain == 0: continue
                _d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean, stop_if_div=False, maxiter=100)
                d_mdl = n.fft.fft(_d) * window + n.fft.fft(info['res'] * area)
                if opts.rmsky: d[:,ch] = (d[:,ch] * window - d_mdl) * w[:,ch]
                else: d[:,ch] = d_mdl
            #dat[bl][pol] = d * w
            dat[bl][pol] = d
            for ti,di,wi in zip(times, dat[bl][pol], window):
                data[bl][pol][ti] = data[bl][pol].get(ti, 0) + di * wi
                wgts[bl][pol][ti] = wgts[bl][pol].get(ti, 0) + wi**2

for filename in args[1:-1]:
    outfile = filename+'X'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print outfile, 'exists.  Skipping...'
        continue
    def mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        #cnt = n.searchsorted(times,t)
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        try:
            d,w = data[bl][pol][t], wgts[bl][pol][t]
            d /= w
        except(KeyError): pass
        return p, d, f
        #return p, d, n.zeros_like(f)
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist='FRINGE RATE FILTER:'+' '.join(sys.argv)+'\n', raw=True)
