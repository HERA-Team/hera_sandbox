#! /usr/bin/env python
"""
Filter in fringe-rate to select (or de-select) fringe rates that correspond to sources fixed
to the celestial sphere.
Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--minfr', type='float', default=6e-5,
    help='Minimum fringe rate (in Hz) to allow.')
o.add_option('--maxfr', type='float', default=6e-5,
    help='Maximum fringe rate (in Hz) to allow.')
o.set_usage('fringe_rate_filter.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

#uv = a.miriad.UV(args[0])
#print "inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly"
#inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#del(uv)

def fng_thresh(inttime, nints, min_fr, max_fr):
    '''Return the (upper,lower) bin indices for the fringe rates indicated.'''
    bin_fr = 1. / (inttime * nints)
    uthr = max_fr / bin_fr
    lthr = min_fr / bin_fr
    uthr, lthr = n.ceil(uthr).astype(n.int), int(n.floor(lthr))
    return uthr,lthr

for filename in args:
    outfile = filename + 'F'
    if os.path.exists(outfile):
        print '    %s exists.  Skipping...' % (outfile)
        continue

    print '    Reading', filename
    times, dat, flg = C.arp.get_dict_of_uv_data([filename], opts.ant, opts.pol, verbose=False)
    inttime = (times[1]-times[0]) / a.ephem.second
    print 'INTTIME (s):', inttime
    #uv = a.miriad.UV(filename)
    ufr,lfr = fng_thresh(inttime, times.size, min_fr=opts.minfr, max_fr=opts.maxfr)

    data,wgts = {}, {}

    def mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        try:
            d_,w = data[bl][pol].pop(t), wgts[bl][pol].pop(t)
            #f_ = n.where(w < .75 * w.max(), 1, 0)
            f_ = f
            d_ /= n.where(f_, n.Inf, w)
        except(KeyError): return p, None, None
        return p, d_, f_


    #window = a.dsp.gen_window(times.shape[0], 'blackman-harris')
    window = 1

    for bl in dat.keys():
        if not data.has_key(bl): data[bl], wgts[bl] = {}, {}
        # Variables: ufr (upper fringe rate), lfr (lowest fringe rate)
        for pol in dat[bl].keys():
            if not data[bl].has_key(pol): data[bl][pol], wgts[bl][pol] = {}, {}
            d = n.where(flg[bl][pol], 0, dat[bl].pop(pol))
            w = n.logical_not(flg[bl].pop(pol)).astype(n.float)
            for ch in range(d.shape[1]):
                _d,_w = n.fft.ifft(d[:,ch]*window), n.fft.ifft(w[:,ch]*window)
                gain = n.abs(_w[0])
                #print ch, gain
                if gain < .01: continue
                area = n.ones(_d.shape, dtype=n.int)
                # XXX would prefer to implement fr cuts as weights rather than bin ranges...
                area[ufr+1:lfr] = 0
                _d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean, stop_if_div=False, maxiter=100)
                _r = info['res'] * area
                #if opts.minfr >= 0:
                #    _d[:lfr+1] = 0
                #    _r[:lfr+1] = 0
                #    if lfr > 0:
                #        _d[-lfr:] = 0
                #        _r[-lfr:] = 0
                w_ch = n.fft.fft(n.fft.ifft(w[:,ch]) * area)#; w_ch = n.where(w_ch > .75, w_ch, 0)
                d_ch = n.fft.fft(_d) * window * w_ch + n.fft.fft(_r)
                d[:,ch] = d_ch
                w[:,ch] = w_ch
            for (ti,di,wi) in zip(times, d, w):
                data[bl][pol][ti] = di
                wgts[bl][pol][ti] = wi
    #filename = files[0]
    #outfile = filename+'F'
    print '    Writing', outfile
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)
    del(uvo) # helps the file be available for reading sooner

