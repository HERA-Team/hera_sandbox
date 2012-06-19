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
o.add_option('--minfr', dest='minfr', type='float', default=6e-5,
    help='Minimum fringe rate (in Hz) to allow.  Anything varying slower than this is considered crosstalk.  A negative value indicates nothing should be considered crosstalk.  Default 6e-5')
o.add_option('--fr_frac', type='float', default=1.,
    help='Fractional width of the fringe-rate filter range to extract.  Default 1.')
o.add_option('--rmsky', action='store_true',
    help='Instead of retaining the data corresponding to the sky, remove it.')
o.set_usage('fringe_rate_filter.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
print "inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly"
inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

def triplets(seq):
    for i in range(1,len(seq)-1): yield seq[i-1:i+2]
    return

def sky_fng_thresh(bl_ew_len, inttime, nints, freq, min_fr=6e-5, neg_fr=-2e-4, max_fr_frac=1.):
    '''For bl_ew_len (the east/west projection) in ns, return the (upper,negative,lower) fringe rate bins 
    that geometrically correspond to the sky.'''
    bin_fr = 1. / (inttime * nints)
    max_bl = bl_ew_len * max_fr_frac
    max_fr = freq * max_bl * 2*n.pi / a.const.sidereal_day
    lthr = min_fr / bin_fr
    nthr = neg_fr / bin_fr
    uthr = max_fr / bin_fr
    uthr, nthr, lthr = n.ceil(uthr).astype(n.int), int(n.floor(nthr)), int(n.floor(lthr))
    return (uthr,nthr,lthr)

def all_sky_fng_thresh(aa, inttime, nints, min_fr=6e-5, neg_fr=-2e-4, max_fr_frac=1.):
    '''Return a dictionary, indexed by baseline, of the (upper,lower) fringe rate
    bins that geometrically correspond to the sky.'''
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl = aa.ij2bl(i,j)
        bl_len = aa.get_baseline(i,j)
        #print i,j, bl_len
        bl_ew_len = n.sqrt(n.dot(bl_len[:2], bl_len[:2]))
        filters[bl] = sky_fng_thresh(bl_ew_len, inttime, nints, aa.get_afreqs(), 
            min_fr=min_fr, neg_fr=neg_fr, max_fr_frac=max_fr_frac)
    return filters

data,wgts = {}, {}

def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    bl = a.miriad.ij2bl(i,j)
    pol = a.miriad.pol2str[uv['pol']]
    try:
        d_,w = data[bl][pol].pop(t), wgts[bl][pol].pop(t)
        f_ = n.where(w < .75 * w.max(), 1, 0)
        d_ /= n.where(f_, n.Inf, w)
    except(KeyError): return p, None, None
    if opts.rmsky: return p, n.where(f_, 0, d-d_), f_
    else: return p, d_, f_

for files in triplets(args):
    print '(%s) %s (%s) ->' % (files[0], files[1], files[2])
    outfiles = [f+'F' for f in files]
    # Don't read triplets for which all output files exist
    if os.path.exists(outfiles[0]) and os.path.exists(outfiles[1]) and os.path.exists(outfiles[2]):
        print '    All output files exist.  Skipping...'
        continue

    # XXX rereading files 3x (for each position in triplet) is inefficient
    # XXX also, this process fails if not all files are the same size...
    print '    Reading files'
    times, dat, flg = C.arp.get_dict_of_uv_data(files, opts.ant, opts.pol, verbose=False)
    uv1 = a.miriad.UV(files[1]); uv2 = a.miriad.UV(files[2])
    t1,t2 = uv1.read()[0][1], uv2.read()[0][1]
    del(uv1); del(uv2)
    times0 = times[n.where(times < t1)]
    times1 = times[n.where(n.logical_and(times < t2, times >= t1))]
    times2 = times[n.where(times >= t2)]
    #print len(times0), len(times1), len(times2)
    max_fr = all_sky_fng_thresh(aa, inttime, times.size, min_fr=opts.minfr, max_fr_frac=opts.fr_frac)

    window = a.dsp.gen_window(times.shape[0], 'blackman-harris')

    for bl in dat.keys():
        if not data.has_key(bl): data[bl], wgts[bl] = {}, {}
        # Variables: ufr (upper fringe rate), nfr (negative fringe rate), lfr (lowest fringe rate for xtalk removal)
        ufr,nfr,lfr = max_fr[bl]
        for pol in dat[bl].keys():
            if not data[bl].has_key(pol): data[bl][pol], wgts[bl][pol] = {}, {}
            d = n.where(flg[bl][pol], 0, dat[bl].pop(pol))
            w = n.logical_not(flg[bl].pop(pol)).astype(n.float)
            for ch in range(d.shape[1]):
                _d,_w = n.fft.ifft(d[:,ch]*window), n.fft.ifft(w[:,ch]*window)
                gain = n.abs(_w[0])
                #print ch, lfr, nfr, ufr[ch], gain
                if gain == 0: continue
                area = n.ones(_d.shape, dtype=n.int)
                # XXX would prefer to implement fr cuts as weights rather than bin ranges...
                area[ufr[ch]+1:nfr] = 0
                _d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean, stop_if_div=False, maxiter=100)
                _r = info['res'] * area
                if opts.minfr >= 0:
                    _d[:lfr+1] = 0
                    _r[:lfr+1] = 0
                    if lfr > 0:
                        _d[-lfr:] = 0
                        _r[-lfr:] = 0
                w_ch = n.fft.fft(n.fft.ifft(w[:,ch]) * area)#; w_ch = n.where(w_ch > .75, w_ch, 0)
                d_ch = n.fft.fft(_d) * window * w_ch + n.fft.fft(_r)
                d[:,ch] = d_ch
                w[:,ch] = w_ch
            for (ti,di,wi,fi) in zip(times, d, window, w):
                # Only process the center file (simpler when decimating)
                if ti in times0 and os.path.exists(outfiles[0]): continue
                if ti in times1 and os.path.exists(outfiles[1]): continue
                if ti in times2 and os.path.exists(outfiles[2]): continue
                data[bl][pol][ti] = data[bl][pol].get(ti, 0) + di * wi * fi
                wgts[bl][pol][ti] = wgts[bl][pol].get(ti, 0) + (wi * fi)**2
    filename = files[0]
    outfile = filename+'F'
    if os.path.exists(outfile):
        print '    %s exists.  Skipping...' % outfile
        continue
    print '    Writing', outfile
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)

for filename in files[1:]:
    outfile = filename+'F'
    if os.path.exists(outfile):
        print '    %s exists.  Skipping...' % outfile
        continue
    print '    Writing', outfile
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)
