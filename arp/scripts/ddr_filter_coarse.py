#! /usr/bin/env python
"""
Generate 3 output files containing components from the original file, organized by smoothness.
Smoothness is determined by the fringe-rate and delay corresponding to the maximum baseline length provided.
..."E" is for smooth-in-freq, smooth-in-time components (nominally, the extracted sky)
..."D" for smooth-in-freq components (delay-filtered)
..."F" for smooth-in-time components (fringe-rate-filtered)
"""

import aipy as a, numpy as n, sys, os, optparse
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--maxbl', type='float', default=300.,
    help='Maximum baseline length (in meters) for purpose of calculating maximum fringe rates. Default 300.')
o.add_option('--lat', type='float', default=-30.,
    help='Latitude of array in degrees.  Default -30.  Used to estimate maximum fringe rates.')
o.add_option('--corrmode', default='j',
    help='Data type ("r" for float32, "j" for shared exponent int16) of dly/fng output files.  Default is "j".')
o.add_option('--invert', action='store_true',
    help='Invert each filter.')
o.add_option('--outputs', default='ddr,dly,fng',
    help='Comma-delimited list of output file types (ddr, dly, fng).  Default ddr,dly,fng (i.e. all of them)')
o.add_option('--no_decimate', action='store_true',
    help='Instead of decimating, match output file size to input file size.')
o.set_usage('ddr_filter_coarse.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

def gen_skypass_delay(bl_len, sdf, nchan):
    bin_dly = 1. / (sdf * nchan)
    uthresh, lthresh = bl_len/bin_dly + 1.5, -bl_len/bin_dly - 0.5
    uthresh, lthresh = int(n.ceil(uthresh)), int(n.floor(lthresh)) 
    return uthresh,lthresh

def sky_fng_thresh(bl_len, inttime, nints, freq, lat):
    '''For bl_ew_len (the east/west projection) in ns, return the (upper,negative,lower) fringe rate bins 
    that geometrically correspond to the sky.'''
    bin_fr = 1. / (inttime * nints)
    max_fr = freq * bl_len * 2*n.pi / a.const.sidereal_day
    neg_fr = max_fr * n.cos(n.pi/2 + n.abs(lat))
    nfng = neg_fr / bin_fr
    ufng = max_fr / bin_fr
    ufng, nfng = n.ceil(ufng).astype(n.int), int(n.floor(nfng))
    return ufng+1, nfng

def triplets(seq):
    for i in range(1,len(seq)-1): yield seq[i-1:i+2]
    return

maxbl = opts.maxbl * 1e2 / a.const.len_ns
uv = a.miriad.UV(args[0])
inttime = uv['inttime']
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
udly,ndly = gen_skypass_delay(maxbl, uv['sdf'], uv['nchan'])
window_dly = a.dsp.gen_window(fqs.shape[0], 'blackman-harris')
_w = n.fft.ifft(window_dly)
_w = n.concatenate([_w[:udly], _w[ndly:]])
window_dly_dec = n.fft.fft(_w)
# Info for creating fqs_dec
nchan_dec = window_dly_dec.size
sdf_dec = uv['sdf'] * uv['nchan'] / nchan_dec
sfreq_dec = uv['sfreq']
# Info for creating times_dec
inttime = uv['inttime']
window_dly.shape = (1,) + window_dly.shape
window_dly_dec.shape = (1,) + window_dly_dec.shape
del(uv)

DECIMATE = not (opts.no_decimate or opts.invert)
opts.outputs = opts.outputs.split(',')
DDR = 'ddr' in opts.outputs
DLY = 'dly' in opts.outputs
FNG = 'fng' in opts.outputs

for files in triplets(args):
    print '(%s) %s (%s) ->' % (files[0], files[1], files[2])
    filename = files[1]
    outfile_dly = filename+'D'
    outfile_fng = filename+'F'
    outfile_ddr = filename+'E'
    if (not DLY or os.path.exists(outfile_dly)) and \
            (not FNG or os.path.exists(outfile_fng)) and \
            (not DDR or os.path.exists(outfile_ddr)):
        print '    All output files exist.  Skipping...'
        continue
    data_ddr, wgts_ddr = {}, {}
    data_dly, wgts_dly = {}, {}
    data_fng, wgts_fng = {}, {}
    match_tdec = {}

    # XXX rereading files 3x (for each position in triplet) is inefficient
    # XXX also, this process fails if not all files are the same size...
    print '    Reading files'
    times, dat, flg = C.arp.get_dict_of_uv_data(files, opts.ant, opts.pol, verbose=False)
    uv1 = a.miriad.UV(files[1]); uv2 = a.miriad.UV(files[2])
    t1,t2 = uv1.read()[0][-1], uv2.read()[0][-1]
    del(uv1); del(uv2)
    times0 = times[n.where(times < t1)]
    times1 = times[n.where(n.logical_and(times < t2, times >= t1))]
    times2 = times[n.where(times >= t2)]
    print len(times0), len(times1), len(times2)
    # Variables: ufng (upper fringe rate), nfng (negative fringe rate)
    ufng,nfng = sky_fng_thresh(maxbl, inttime, len(times), fqs.max(), opts.lat*a.img.deg2rad)
    # Make fringe filter width divisible by 3 so that decimation pattern is periodic across files
    if (ufng-nfng) % 3 == 1: ufng,nfng = ufng+1,nfng-1
    elif (ufng-nfng) % 3 == 2: ufng = ufng + 1
    assert((ufng-nfng) % 3 == 0)

    window_fng = a.dsp.gen_window(times.shape[0], 'blackman-harris')
    _w = n.fft.ifft(window_fng)
    _w = n.concatenate([_w[:ufng], _w[nfng:]])
    window_fng_dec = n.fft.fft(_w)
    window_fng.shape = window_fng.shape + (1,)
    window_fng_dec.shape = window_fng_dec.shape + (1,)
    window = window_fng * window_dly
    # Create times_dec
    start_jd = times[0]
    nints = times.size
    nints_dec = window_fng_dec.size
    delta_jd = n.average(times[1:] - times[:-1])
    if DECIMATE:
        times_dec = start_jd + n.arange(nints_dec,dtype=n.float) * delta_jd * nints / nints_dec
        inttime_dec = inttime * nints / nints_dec
    else:
        times_dec = times
        inttime_dec = inttime
    # Match each dec time to a single original time for use in mfunc later
    #print nints, nints_dec
    for cnt,tdec in enumerate(times_dec):
        if cnt < times_dec.size / 3 or cnt >= 2 * times_dec.size / 3: continue
        t = times[n.argmin(n.abs(times-tdec))]
        match_tdec[t] = tdec
    # Create area for allowable clean components
    area_dly = n.ones(times.shape + fqs.shape, dtype=n.int)
    area_dly[:,udly:ndly] = 0
    area_fng = n.ones(times.shape + fqs.shape, dtype=n.int)
    area_fng[ufng:nfng,:] = 0
    area = area_dly * area_fng
    #print (ufng,nfng,udly,ndly)
    print '    Processing data'
    for bl in dat:
        if not data_ddr.has_key(bl):
            # XXX the memory footprint of all these copies of data is high.  need to reduce
            data_ddr[bl],wgts_ddr[bl] = {}, {}
            data_fng[bl],wgts_fng[bl] = {}, {}
            data_dly[bl],wgts_dly[bl] = {}, {}
        for pol in dat[bl].keys():
            if not data_ddr[bl].has_key(pol):
                data_ddr[bl][pol],wgts_ddr[bl][pol] = {}, {}
                data_fng[bl][pol],wgts_fng[bl][pol] = {}, {}
                data_dly[bl][pol],wgts_dly[bl][pol] = {}, {}
            d = n.where(flg[bl][pol], 0, dat[bl][pol])
            del(dat[bl][pol]); del(flg[bl][pol]) # We're done with the raw data
            w = n.where(n.abs(d) == 0, 0., 1.)
            _d,_w = n.fft.ifft2(d*window), n.fft.ifft2(w*window)
            gain = n.abs(_w[0,0])
            if gain == 0: continue
            _d,info = a.deconv.clean(_d,_w, area=area, tol=opts.clean, stop_if_div=False, maxiter=100)
            _r = info['res']
            if DECIMATE:
                _f = n.fft.ifft2(w)
                # Filter and decimate for dly, fng, and ddr datasets
                data_L,wgts_L,t_L,d_L,w_L,f_L = [],[],[],[],[],[]
                if DDR:
                    _d_ddr = n.concatenate([_d[:ufng], _d[nfng:]])
                    _d_ddr = n.concatenate([_d_ddr[:,:udly], _d_ddr[:,ndly:]], axis=1)
                    _r_ddr = n.concatenate([_r[:ufng], _r[nfng:]])
                    _r_ddr = n.concatenate([_r_ddr[:,:udly], _r_ddr[:,ndly:]], axis=1)
                    _f_ddr = n.concatenate([_f[:ufng], _f[nfng:]])
                    _f_ddr = n.concatenate([_f_ddr[:,:udly], _f_ddr[:,ndly:]], axis=1)
                _r *= n.logical_not(area) # Remove the part that went in the DDR file from all the rest
                if FNG:
                    _r_fng = n.concatenate([_r[:ufng], _r[nfng:]])
                    _f_fng = n.concatenate([_f[:ufng], _f[nfng:]])
                    d_fng = n.fft.fft2(_r_fng)
                    w_fng = window_fng_dec * window_dly
                    f_fng = n.fft.fft2(_f_fng); f_fng = n.where(f_fng > .75, f_fng, 0)
                    t_fng = times_dec
                    data_L.append(data_fng); wgts_L.append(wgts_fng)
                    t_L.append(t_fng); d_L.append(d_fng); w_L.append(w_fng); f_L.append(f_fng)
                if DLY:
                    _r_dly = n.concatenate([_r[:,:udly], _r[:,ndly:]], axis=1)
                    _f_dly = n.concatenate([_f[:,:udly], _f[:,ndly:]], axis=1)
                    d_dly = n.fft.fft2(_r_dly)
                    w_dly = window_fng * window_dly_dec
                    f_dly = n.fft.fft2(_f_dly); f_dly = n.where(f_dly > .75, f_dly, 0)
                    t_dly = times
                    data_L.append(data_dly); wgts_L.append(wgts_dly)
                    t_L.append(t_dly); d_L.append(d_dly); w_L.append(w_dly); f_L.append(f_dly)
                if DDR:
                    #d_ddr = n.fft.fft2(_d_ddr) * (window_fng_dec * window_dly_dec) + n.fft.fft2(_r_ddr)
                    w_ddr = window_fng_dec * window_dly_dec
                    f_ddr = n.fft.fft2(_f_ddr); f_ddr = n.where(f_ddr > .75, f_ddr, 0)
                    d_ddr = n.fft.fft2(_d_ddr) * w_ddr * f_ddr + n.fft.fft2(_r_ddr)
                    t_ddr = times_dec
                    data_L.append(data_ddr); wgts_L.append(wgts_ddr)
                    t_L.append(t_ddr); d_L.append(d_ddr); w_L.append(w_ddr); f_L.append(f_ddr)
                # Whichever it was (decimate or not) collect the data into dictionaries
                for data_dat, wgts_dat, t_dat, d_dat, w_dat, f_dat in zip(data_L, wgts_L, t_L, d_L, w_L, f_L):
                    for cnt,(ti,di,wi,fi) in enumerate(zip(t_dat, d_dat, w_dat, f_dat)):
                        # Only process the center file (simpler when decimating)
                        if cnt < t_dat.size / 3 or cnt >= 2 * t_dat.size / 3: continue
                        data_dat[bl][pol][ti] = data_dat[bl][pol].get(ti, 0) + di * wi * fi
                        wgts_dat[bl][pol][ti] = wgts_dat[bl][pol].get(ti, 0) + (wi * fi)**2
            else: # Don't decimate
                if DLY:
                    d_dly = n.fft.fft2(_r * area_dly * n.logical_not(area))
                    f_dly = n.fft.fft2(n.fft.ifft2(w) * area_dly); f_dly = n.where(f_dly > .75, f_dly, 0)
                    for cnt,(ti,di,wi,fi) in enumerate(zip(times, d_dly, window, f_dly)):
                        # Only process the center file (simpler when decimating)
                        if cnt < times.size / 3 or cnt >= 2 * times.size / 3: continue
                        data_dly[bl][pol][ti] = data_dly[bl][pol].get(ti, 0) + di * wi * fi
                        wgts_dly[bl][pol][ti] = wgts_dly[bl][pol].get(ti, 0) + (wi * fi)**2
                if FNG:
                    d_fng = n.fft.fft2(_r * area_fng * n.logical_not(area))
                    f_fng = n.fft.fft2(n.fft.ifft2(w) * area_fng); f_fng = n.where(f_fng > .75, f_fng, 0)
                    for cnt,(ti,di,wi,fi) in enumerate(zip(times, d_fng, window, f_fng)):
                        # Only process the center file (simpler when decimating)
                        if cnt < times.size / 3 or cnt >= 2 * times.size / 3: continue
                        data_fng[bl][pol][ti] = data_fng[bl][pol].get(ti, 0) + di * wi * fi
                        wgts_fng[bl][pol][ti] = wgts_fng[bl][pol].get(ti, 0) + (wi * fi)**2
                if DDR:
                    #d_ddr = n.fft.fft2(_d) * window + n.fft.fft2(_r * area)
                    #w_ddr = n.fft.fft2(_w * area)
                    f_ddr = n.fft.fft2(n.fft.ifft2(w) * area); f_ddr = n.where(f_ddr > .75, f_ddr, 0)
                    d_ddr = n.fft.fft2(_d) * window * f_ddr + n.fft.fft2(_r * area)
                    for cnt,(ti,di,wi,fi) in enumerate(zip(times, d_ddr, window, f_ddr)):
                        # Only process the center file (simpler when decimating)
                        if cnt < times.size / 3 or cnt >= 2 * times.size / 3: continue
                        data_ddr[bl][pol][ti] = data_ddr[bl][pol].get(ti, 0) + di * wi * fi
                        wgts_ddr[bl][pol][ti] = wgts_ddr[bl][pol].get(ti, 0) + (wi * fi)**2
    def mfunc_dly(uv, p, d, f):
        uvw,t,(i,j) = p
        p = uvw,t,(i,j)
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        try:
            d_,w = data_dly[bl][pol][t], wgts_dly[bl][pol][t]
            d_ /= n.where(w > 0, w, 1)
            f_ = n.where(w == 0, 1, 0)
        except(KeyError): return p, None, None
        if opts.invert: return p, n.where(f_, 0, d-d_), f_
        else: return p, d_, f_
    def mfunc_fng(uv, p, d, f):
        uvw,t,(i,j) = p
        try: t = match_tdec[t]
        except(KeyError): return p, None, None
        p = uvw,t,(i,j)
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        try:
            d_,w = data_fng[bl][pol][t], wgts_fng[bl][pol][t]
            d_ /= n.where(w > 0, w, 1)
            f_ = n.where(w == 0, 1, 0)
        except(KeyError): return p, None, None
        if opts.invert: return p, n.where(f_, 0, d-d_), f_
        else: return p, d_, f_
    def mfunc_ddr(uv, p, d, f):
        uvw,t,(i,j) = p
        try: t = match_tdec[t]
        except(KeyError): return p, None, None
        p = uvw,t,(i,j)
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        try:
            d_,w = data_ddr[bl][pol][t], wgts_ddr[bl][pol][t]
            d_ /= n.where(w > 0, w, 1)
            f_ = n.where(w == 0, 1, 0)
        except(KeyError): return p, None, None
        if opts.invert: return p, n.where(f_, 0, d-d_), f_
        else: return p, d_, f_
    uvi = a.miriad.UV(filename)

    if DLY:
        if not os.path.exists(outfile_dly):
            print '   ', outfile_dly
            # For noise-like files (delay, fringe filter), usually use data type 'j' for enhanced compression.
            uvo = a.miriad.UV(outfile_dly, corrmode=opts.corrmode, status='new')
            if DECIMATE: uvo.init_from_uv(uvi, override={'sdf':sdf_dec, 'sfreq':sfreq_dec, 'nchan':nchan_dec})
            else: uvo.init_from_uv(uvi)
            uvo.pipe(uvi, mfunc=mfunc_dly, append2hist='DDR FILTER, DLY:'+' '.join(sys.argv)+'\n', raw=True)
            del(uvo)
        else: print '   ', outfile_dly, 'exists.  Skipping...'
        uvi.rewind()

    if FNG:
        if not os.path.exists(outfile_fng):
            print '   ', outfile_fng
            # For noise-like files (delay, fringe filter), usually use data type 'j' for enhanced compression.
            uvo = a.miriad.UV(outfile_fng, corrmode=opts.corrmode, status='new')
            if DECIMATE: uvo.init_from_uv(uvi, override={'sdf':sdf_dec, 'sfreq':sfreq_dec, 'inttime':inttime_dec})
            else: uvo.init_from_uv(uvi)
            uvo.pipe(uvi, mfunc=mfunc_fng, append2hist='DDR FILTER, FNG:'+' '.join(sys.argv)+'\n', raw=True)
            del(uvo)
        else: print '   ', outfile_fng, 'exists.  Skipping...'
        uvi.rewind()

    if DDR:
        if not os.path.exists(outfile_ddr):
            print '   ', outfile_ddr
            uvo = a.miriad.UV(outfile_ddr, status='new')
            if DECIMATE: uvo.init_from_uv(uvi, override={'sdf':sdf_dec, 'sfreq':sfreq_dec, 'nchan':nchan_dec, 'inttime':inttime_dec})
            else: uvo.init_from_uv(uvi)
            uvo.pipe(uvi, mfunc=mfunc_ddr, append2hist='DDR FILTER, DDR:'+' '.join(sys.argv)+'\n', raw=True)
            del(uvo)
        else: print '   ', outfile_ddr, 'exists.  Skipping...'
