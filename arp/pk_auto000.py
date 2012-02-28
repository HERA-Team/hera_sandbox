#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
o.set_usage('uv2pk.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True)
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
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
MIN_CH,MAX_CH,SUBBAND = 150,900,50
#MIN_CH,MAX_CH,SUBBAND = 300,1800,50
NTAPS = 3

SCALE = 30e3
FoV = 1.
X2Y = 730.

DECONV1WIDTH = 5
PLOT = True

def dk_deta(z):
    '''(1/Mpc)/(1/MHz)'''
    return 5.1 * 0.073 * ((1+z) / 10.)**-.5

_sum,_wgt = {},{}
times = []
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    src = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
        if not _sum.has_key(i): _sum[i],_wgt[i] = {},{}
        d = d[MIN_CH:MAX_CH] * SCALE
        #if (len(times) + i) % 2 == 0: d *= -1
        #if (len(times)) % 4 < 2: d *= -1
        #if (len(times)) >= 1344/2: d *= -1
        val = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        d *= val
        if True:  # This somewhat reduces src sidelobes
            #sys.stdout.write('.'); sys.stdout.flush()
            gain = n.sqrt(n.average(val**2))
            if gain == 0: continue
            if False: mdl = (freqs/.150)**-2.5
            elif True:
                poly = n.polyfit(n.log10(freqs[MIN_CH:MAX_CH].compress(val)), n.log10(n.abs(d.compress(val))), deg=3)
                mdl = 10**(n.polyval(poly, n.log10(freqs)))
            else: mdl = n.ones_like(freqs)
            ker = n.fft.ifft(val*mdl[MIN_CH:MAX_CH])
            _d = n.fft.ifft(d)
            __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            #__d[1:] = 0; _d -= n.fft.ifft(n.fft.fft(__d) * val)
            if DECONV1WIDTH == 0: __d[1:] = 0
            else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
            if PLOT:
                if False:
                    p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(d))
                    p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(n.fft.fft(__d) * n.fft.fft(ker)))
                    p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(d - n.fft.fft(__d) * n.fft.fft(ker)))
                    p.ylim(1e2,1e6)
                else:
                    p.plot(freqs[MIN_CH:MAX_CH], d)
                    p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(__d) * n.fft.fft(ker))
                    #p.plot(freqs[MIN_CH:MAX_CH], d - n.fft.fft(__d) * n.fft.fft(ker))
            _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(ker))
            if True:
                ker = n.fft.ifft(val)
                __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
                if DECONV1WIDTH == 0: __d[1:] = 0
                else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
                if PLOT:
                    #p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(__d) * n.fft.fft(ker))
                    p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(_d) - n.fft.fft(__d) * n.fft.fft(ker))
                    _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(ker))
            d = n.fft.fft(_d) * val
        w = val
        for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
            ch0 = b*SUBBAND/2
            ch1 = ch0 + NTAPS*SUBBAND
            if ch1 > MAX_CH-MIN_CH: continue
            #if not _sum.has_key(b): _sum[b],_wgt[b] = {},{}
            _d = pfb(d[ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
            _w = pfb(w[ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
            _sum[i][b] = _sum[i].get(b,0) + _d
            _wgt[i][b] = _wgt[i].get(b,0) + _w
    if PLOT: p.show()

data,gain = {},{}
__sum,__wgt,__dec,__cln = {},{},{},{}
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    #_freqs = freqs[ch0+int((NTAPS-1)/2.):ch0+int((NTAPS-1)/2.)+SUBBAND]
    #B = (_freqs[-1] - _freqs[0]) * 1e9
    #z = n.average(1.420405 / _freqs - 1)
    #eta = n.arange(0, 1./(_freqs[1]-_freqs[0]), 1./(_freqs[-1]-_freqs[0]))
    #eta = n.where(eta > .5/(_freqs[1]-_freqs[0]), eta - (1./(_freqs[1]-_freqs[0])), eta)
    #k_pl = (eta * 1e-3 * dk_deta(z))[:d.size/2]
    #scalar = X2Y * k_pl**3 / (2*n.pi**2) * FoV * B / 8
    for i in _sum:
        g = n.abs(_wgt[i][b][0])
        if True:
            __d, info = a.deconv.clean(_sum[i][b], _wgt[i][b], tol=opts.clean)
            __d[1:] = 0
            #p.plot(n.fft.fft(_sum[i][b])/g)
            _sum[i][b] -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(_wgt[i][b]))
            #p.plot(n.fft.fft(_sum[i][b])/g)
            #p.show()
        d = _sum[i][b][:_sum[i][b].size/2]
        print i, g
        for _d,_g in zip(data.get(b,[]), gain.get(b,[])):
            __sum[b] = __sum.get(b,0) + d * n.conj(_d)
            __wgt[b] = __wgt.get(b,0) + g * _g
            #__wgt.append(n.abs(_wgt[b][:d.size/2]/gain)**2)
            #__dec.append(n.abs(info['res'][:d.size/2]/gain)**2 * scalar)
            #__cln.append(n.abs(d[:d.size/2])**2)
        data[b] = data.get(b,[]) + [d]
        gain[b] = gain.get(b,[]) + [g]

dat = []
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    _freqs = freqs[ch0+int((NTAPS-1)/2.):ch0+int((NTAPS-1)/2.)+SUBBAND]
    B = (_freqs[-1] - _freqs[0]) * 1e9
    z = n.average(1.420405 / _freqs - 1)
    eta = n.arange(0, 1./(_freqs[1]-_freqs[0]), 1./(_freqs[-1]-_freqs[0]))
    eta = n.where(eta > .5/(_freqs[1]-_freqs[0]), eta - (1./(_freqs[1]-_freqs[0])), eta)
    k_pl = (eta * 1e-3 * dk_deta(z))[:__sum[b].size]
    #scalar = X2Y * k_pl**3 / (2*n.pi**2) * FoV * B / 8
    scalar = 1
    dat.append(n.abs(__sum[b] / __wgt[b]) * scalar)

#p.subplot(221)
dat1 = n.log10(dat)
#mx,drng = 10,7
#mx,drng = 7,5
mx,drng = 4,3
#mx,drng = dat1.max(), dat1.max()-max(dat1.min(),2)
p.imshow(dat1, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
p.colorbar(shrink=.5)
#p.subplot(222)
#dat2 = n.log10(__wgt)
#p.imshow(dat2, vmax=dat2.max(), vmin=dat2.max()-8, aspect='auto', interpolation='nearest')
#p.colorbar(shrink=.5)
#p.subplot(223)
#dat3 = n.log10(__dec)
#p.imshow(dat3, vmax=mx, vmin=mx-4, aspect='auto', interpolation='nearest')
#p.colorbar(shrink=.5)
#p.subplot(224)
#dat4 = n.log10(__cln)
#p.imshow(dat4, vmax=dat4.max(), vmin=dat4.max()-8, aspect='auto', interpolation='nearest')
#p.colorbar(shrink=.5)
print len(times)
p.show()


