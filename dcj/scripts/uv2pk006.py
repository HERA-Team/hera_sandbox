#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('uv2pk.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

pm = {'L': None, 'taps':None, 'fwidth':None, 'window':None, 'window_name':None, 'sinx_x':None, 'window_sinx_x':None}

def __set_pm__(L, window, taps, fwidth):
    global pm
    if pm['L'] == L and pm['taps'] == taps and pm['fwidth'] == fwidth:
        if type(window) == str and pm['window_name'] == window: return
        elif window is pm['window']: return
    else:
        pm['L'] = L
        pm['taps'] = taps
        pm['fwidth'] = fwidth
        def sinx_x(x):
            t = n.pi * taps * fwidth * (x/float(L) - .5)
            v = n.where(t != 0, t, 1)
            return n.where(t != 0, n.sin(v) / v, 1)
        pm['sinx_x'] = n.fromfunction(sinx_x, (L,))
    if type(window) == str:
        wf = {}
        wf['hamming'] = lambda x: .54 + .46 * cos(2*n.pi*x/L - n.pi)
        wf['hanning'] = lambda x: .5 + .5 * n.cos(2*n.pi*x/(L+1) - n.pi)
        wf['none'] = lambda x: 1
        pm['window'] = n.fromfunction(wf[window], (L,))
        pm['window_name'] = window
    else:
        pm['window'] = window
        pm['window_name'] = None
    pm['window_sinx_x'] = pm['window'] * pm['sinx_x']

def __pfb_fir__(data, window='hamming', taps=8, fwidth=1):
    L = data.shape[-1]
    __set_pm__(L, window, taps, fwidth)
    d = data * pm['window_sinx_x']
    try: d.shape = d.shape[:-1] + (taps, L/taps)
    except: raise ValueError("More taps than samples")
    return n.sum(d, axis=len(d.shape) - 2)

def pfb(data, window='hamming', taps=8, fwidth=1, fft=n.fft.fft):
    """Perform PFB on last dimension of 'data' for multi-dimensional arrays.
    'window' may be a name (e.g. 'hamming') or an array with length of the
    last dimension of 'data'.  'taps' is the number of PFB taps to use.  The
    number of channels out of the PFB will be length out of the last 
    dimension divided by the number of taps. 'fwidth' scales the width of 
    each channel bandpass (to create overlapping filters, for example)."""
    return fft(__pfb_fir__(data, window, taps, fwidth))

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])

#MIN_CH,MAX_CH,SUBBAND = 60,804,100
MIN_CH,MAX_CH,SUBBAND = 60,804,50
NTAPS = 2
aa.select_chans(n.arange(MIN_CH, MAX_CH))
UV_RES = 1.5
LST_RES = 2*n.pi/36

def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(LST_RES)
    return u,v, lst

freqs = aa.get_afreqs()
lstbins = n.arange(0,2*n.pi, LST_RES)

for filename in args:
    print 'Reading', filename
    exists = False
    for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
        ch0 = b*SUBBAND/2
        ch1 = ch0 + SUBBAND
        fq = n.average(freqs[ch0:ch1])
        ofile = '%s__%5.3f.npz' % (filename, fq)
        if os.path.exists(ofile):
            print ofile, 'exists.  Skipping...'
            exists = True
            break
    if exists: continue
    _sum,_wgt = {},{}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    times = []
    src = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            aa.set_jultime(t)
            # need to bin in lst to accumulate coherently...
            lst = lstbins[n.argmin(n.abs(lstbins - aa.sidereal_time()))]
            src = a.phs.RadioFixedBody(lst, aa.lat)
            src.compute(aa)
            times.append(t)
            bm_wgt = aa[0].bm_response(src.get_crds('top'), pol='y')**2
            bm_wgt = bm_wgt.flatten()
        bl = a.miriad.ij2bl(i,j)
        u,v,w = aa.gen_uvw(i,j, src=src)
        u,v = u.flatten(), v.flatten()
        d = d[MIN_CH:MAX_CH]
        val = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        # mask out overly short baselines adversely affected by xtalk removal
        val *= n.where(n.abs(u) < 12.5, 0, 1)
        d = aa.phs2src(d, src, i, j)
        d /= aa.passband(i,j)
        d *= val
        if n.average(val) < .5: continue
        if True:  # This somewhat reduces src sidelobes
            #sys.stdout.write('.'); sys.stdout.flush()
            gain = n.sqrt(n.average(val**2))
            ker = n.fft.ifft(val)
            _d = n.fft.ifft(d)
            _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            if True: _d = info['res']
            else: _d += info['res'] / gain
            d = n.fft.fft(_d) * val
        # Weight by square of primary beam response in the direction we're projecting for.
        d *= bm_wgt
        w = bm_wgt**2 * val
        for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
            ch0 = b*SUBBAND/2
            ch1 = ch0 + NTAPS*SUBBAND
            if ch1 > u.size: continue
            _u = n.average(u[ch0:ch1]) # Doesn't matter that ch1-ch0 > SUBBAND
            _v = n.average(v[ch0:ch1]) # Doesn't matter that ch1-ch0 > SUBBAND
            bin = uv2bin(_u,_v,lst,uv_res=UV_RES,lst_res=LST_RES)
            if not _sum.has_key(b): _sum[b],_wgt[b] = {},{}
            _d = pfb(d[ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
            _w = pfb(w[ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
            _sum[b][bin] = _sum[b].get(bin,0) + _d
            _wgt[b][bin] = _wgt[b].get(bin,0) + _w

    for b in _sum:
        ch0 = b*SUBBAND/2
        ch1 = ch0 + NTAPS*SUBBAND
        fq = n.average(freqs[ch0:ch1]) # Doesn't matter that ch1-ch0 > SUBBAND
        ofile = '%s__%5.3f.npz' % (filename, fq)
        d = {'freqs':freqs[ch0:ch1], 'uvres':UV_RES, 'lstres':LST_RES}
        for bin in _sum[b]:
            d['sum_%d' % (bin)] = _sum[b][bin]
            d['wgt_%d' % (bin)] = _wgt[b][bin]
        print '   Writing', ofile
        n.savez(ofile, **d)

