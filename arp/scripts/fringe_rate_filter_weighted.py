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
o.set_usage('fringe_rate_filter.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
print "inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly"
inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly
print "inttime =", inttime
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
afreqs = aa.get_afreqs()
del(uv)

#def fng_wgts(fng_bin, fq, bl_ew_len, cen164=.000915, wid164=.000160): # old one based off of fng_bm**2, which is incorrect
def fng_wgts(fng_bin, fq, bl_ew_len, cen164=.000908, wid164=.000253):
    cen = cen164 * (fq * bl_ew_len / 16.43) # Calibrated to .164 GHz, 30m baseline
    wid = wid164 * (fq * bl_ew_len / 16.43) # Calibrated to .164 GHz, 30m baseline
    max_fr = fq * bl_ew_len * 2*n.pi / a.const.sidereal_day
    min_fr = -0.5 * fq * bl_ew_len * 2*n.pi / a.const.sidereal_day
    #return n.exp(-(fng_bin-cen)**2/(2*wid**2))
    return n.exp(-(fng_bin-cen)**2/(2*wid**2)) * n.where(fng_bin > max_fr, n.exp(-(fng_bin-max_fr)**2/(2*.00005**2)), 1)
    #return n.where(fng_bin > max_fr, 0, n.where(fng_bin < min_fr, 0, 1)) # XXX shouldn't change sig amplitude
    #return n.where(fng_bin > 1.5*max_fr, 0, n.where(fng_bin < 1.5*min_fr, 0, 1)) # XXX shouldn't change sig amplitude

def triplets(seq):
    for i in range(1,len(seq)-1): yield seq[i-1:i+2]
    return

def quintuplets(seq):
    for i in range(2,len(seq)-2): yield seq[i-2:i+3]
    return

def sky_fng_thresh(bl_ew_len, inttime, nints, freq, min_fr_frac=-.3, xtalk=-1, max_fr_frac=1.):
    '''For bl_ew_len (the east/west projection) in ns, return the (upper,negative,lower) fringe rate bins 
    that geometrically correspond to the sky.'''
    bin_fr = 1. / (inttime * nints)
    max_bl = bl_ew_len * max_fr_frac
    min_bl = bl_ew_len * min_fr_frac
    max_fr = freq * max_bl * 2*n.pi / a.const.sidereal_day
    min_fr = freq * min_bl * 2*n.pi / a.const.sidereal_day
    lthr = xtalk / bin_fr
    uthr = max_fr / bin_fr
    nthr = min_fr / bin_fr
    uthr, nthr, lthr = n.ceil(uthr).astype(n.int), n.ceil(nthr).astype(n.int), int(n.floor(lthr))
    return (uthr,nthr,lthr)

def all_sky_fng_thresh(aa, inttime, nints, min_fr_frac=-.3, xtalk=-1, max_fr_frac=1.):
    '''Return a dictionary, indexed by baseline, of the (upper,lower) fringe rate
    bins that geometrically correspond to the sky.'''
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl = aa.ij2bl(i,j)
        bl_len = aa.get_baseline(i,j)
        #print i,j, bl_len
        #bl_ew_len = n.sqrt(n.dot(bl_len[:2], bl_len[:2]))
        bl_ew_len = bl_len[0] # XXX this assumes the ey component of the bl is ~0
        filters[bl] = sky_fng_thresh(bl_ew_len, inttime, nints, aa.get_afreqs(), 
            min_fr_frac=min_fr_frac, xtalk=xtalk, max_fr_frac=max_fr_frac)
        filters[bl] = (filters[bl], bl_ew_len)
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
    #if opts.rmsky: return p, n.where(f_, 0, d-d_), f_
    else: return p, d_, f_

#for files in triplets(args):
for files in quintuplets(args):
    #print '(%s) %s (%s) ->' % (files[0], files[1], files[2])
    print '((%s)) (%s) %s (%s) ((%s))->' % (files[0], files[1], files[2], files[3], files[4])
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
    fng_bins = n.fft.fftfreq(times.size, inttime)
    max_fr = all_sky_fng_thresh(aa, inttime, times.size, 
            min_fr_frac=-0.75, xtalk=-1, max_fr_frac=1.5)

    window = a.dsp.gen_window(times.shape[0], 'blackman-harris')

    for bl in dat.keys():
        if not data.has_key(bl): data[bl], wgts[bl] = {}, {}
        # Variables: ufr (upper fringe rate), nfr (negative fringe rate), lfr (lowest fringe rate for xtalk removal)
        ufr,nfr,lfr = max_fr[bl][0]
        bl_ew_len = max_fr[bl][1]
        for pol in dat[bl].keys():
            if not data[bl].has_key(pol): data[bl][pol], wgts[bl][pol] = {}, {}
            d = n.where(flg[bl][pol], 0, dat[bl].pop(pol))
            w = n.logical_not(flg[bl].pop(pol)).astype(n.float)
            for ch in range(d.shape[1]):
                #ufr[ch] = ufr[d.shape[1]-1] # XXX
                #nfr[ch] = nfr[d.shape[1]-1] # XXX
                _d,_w = n.fft.ifft(d[:,ch]*window), n.fft.ifft(w[:,ch]*window)
                gain = n.abs(_w[0])
                #print ch, lfr, nfr, ufr[ch], gain
                if gain == 0: continue
                area = n.ones(_d.shape, dtype=n.int)
                warea = n.ones_like(area)
                #print ch, ufr[ch], nfr[ch]
                if ufr[ch] >= 0:
                    if nfr[ch] < 0:
                        area[ufr[ch]+1:nfr[ch]] = 0
                        warea[ufr[ch]+1:nfr[ch]] = 0
                    else:
                        area[:nfr[ch]] = 0
                        area[ufr[ch]+1:] = 0
                        warea[ufr[ch]+1-nfr[ch]:] = 0 # Same width but includes DC term for flagging purposes
                else:
                    if nfr[ch] < 0:
                        area[nfr[ch]:] = 0
                        area[:ufr[ch]] = 0
                        warea[-nfr[ch]-ufr[ch]:] = 0 # Same width but includes DC term for flagging purposes
                    else:
                        area[nfr[ch]+1:ufr[ch]] = 0
                        warea[nfr[ch]+1:ufr[ch]] = 0
                if ch == 130 and bl in [a.miriad.ij2bl(0,16), a.miriad.ij2bl(8,16)]: print ufr[ch], nfr[ch]
                #if bl == a.miriad.ij2bl(0,16): print ch, area.size, area.sum()
                #_d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean, stop_if_div=False, maxiter=100)
                _d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean)
                # do beam weighting of fringe-rate bins
                _d_fwgts = fng_wgts(fng_bins, afreqs[ch], bl_ew_len)
                #_w_fwgts = fng_wgts(fng_bins, afreqs[ch], bl_ew_len, cen164=0) # same width, centered at 0 for DC term
                #_d_fwgts = fng_wgts(fng_bins, .150, bl_ew_len)
                #_w_fwgts = fng_wgts(fng_bins, .150, bl_ew_len, cen164=0) # same width, centered at 0 for DC term
                #win_win = n.concatenate([n.zeros(times0.size), n.ones(times1.size), n.zeros(times2.size)])
                if ch == 100: print 43. / (_d_fwgts.sum() / _d_fwgts.size)
                d_fir_coeffs = n.fft.fft(_d_fwgts)
                #d_fir_coeffs = n.fft.ifft(_d_fwgts)
                
                if False:
                    SZ = d_fir_coeffs.size
                    print SZ
                    d_fir_coeffs = n.concatenate([d_fir_coeffs[-SZ/2+1:], d_fir_coeffs[:-SZ/2+1]])
                    #d_fir_coeffs = d_fir_coeffs[:SZ] * a.dsp.gen_window(SZ, 'blackman-harris')
                else:
                    d_fir_coeffs = n.fft.fftshift(d_fir_coeffs)
                    #print d_fir_coeffs.sum(), n.sum(n.abs(d_fir_coeffs)**2), d_fir_coeffs.size
                    #d_fir_coeffs /= d_fir_coeffs.sum()
                    d_fir_coeffs /= d_fir_coeffs.size
                w_fir_coeffs = n.abs(d_fir_coeffs)
                if False and ch == 100:
                    print d_fir_coeffs.sum(), n.abs(d_fir_coeffs).sum(), d_fir_coeffs.max()
                    import pylab
                    pylab.subplot(211)
                    pylab.plot(d_fir_coeffs.real)
                    pylab.plot(d_fir_coeffs.imag)
                    pylab.plot(w_fir_coeffs)
                    pylab.subplot(212); pylab.plot(_d_fwgts)
                    pylab.show()
                #print n.fft.fft(_d_fwgts)
                #import pylab; pylab.plot(n.fft.fft(_d_fwgts)); pylab.show()
                #_d_fwgts = n.ones_like(_d_fwgts) # XXX
                #_w_fwgts = n.ones_like(_w_fwgts) # XXX
                _r = info['res'] * area
                #if opts.xtalk >= 0:
                if False:
                    _d[:lfr+1] = 0
                    _r[:lfr+1] = 0
                    if lfr > 0:
                        _d[-lfr:] = 0
                        _r[-lfr:] = 0
                #w_ch = n.fft.fft(n.fft.ifft(w[:,ch]) * area)#; w_ch = n.where(w_ch > .75, w_ch, 0)
                # use warea instead of area to include DC term
                #w_ch = n.fft.fft(n.fft.ifft(w[:,ch]) * warea * _w_fwgts)
                #d_ch = n.fft.fft(_d * _d_fwgts) * window * w_ch + n.fft.fft(_r * _d_fwgts)
                w_ch = n.convolve(w[:,ch]*window, w_fir_coeffs, 'same') # XXX is this a problem?  should it be ones?
                #w_ch = w_ch[times1.size:-times1.size] # XXX
                d_ch = n.convolve(d[:,ch]*window, d_fir_coeffs, 'same')
                #d_ch = n.conj(n.convolve(n.conj(d[:,ch])*window, d_fir_coeffs, 'same')) # XXX hack for reverse-oriented baselines
                #d_ch = d_ch[times1.size:-times1.size] # XXX
                if False and ch == 100:
                    print times0.size, times1.size, times2.size
                    import pylab as p
                    #p.plot(d[:,ch], 'k')
                    p.plot(d[:,ch]*window, 'k:')
                    #p.plot(d_ch/w_ch.max(), 'g')
                    #p.plot(w_ch/w_ch.max(), 'b')
                    #p.plot(n.fft.fft(_d_fwgts), 'k')
                    #p.plot(fir_win, 'r')
                    #p.plot(d_ch/w_ch/n.sqrt(n.average(window**2)), 'r')
                    p.plot(d_ch, 'r:')
                    p.show()
                d[:,ch] = d_ch
                #w[:,ch] = w_ch # XXX
            for (ti,di,wi,fi) in zip(times, d, window, w):
            #for (ti,di,wi,fi) in zip(times, d, n.ones_like(window), w):
                #if ti in times0 or ti in times2: continue# XXX Only process the center file (simpler when decimating)
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
    del(uvo) # helps the file be available for reading sooner

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
