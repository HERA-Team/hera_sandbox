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

#DECONV1WIDTH = 5
DECONV1WIDTH = 10
PLOT1 = False
PLOT2 = True

def dk_deta(z):
    '''(1/Mpc)/(1/MHz)'''
    return 5.1 * 0.073 * ((1+z) / 10.)**-.5

dat,wgt = {},{}
times = []
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    src = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
        d = d[MIN_CH:MAX_CH] * SCALE
        #if (len(times) + i) % 2 == 0: d *= -1
        #if (len(times)) % 4 < 2: d *= -1
        #if (len(times)) >= 1344/2: d *= -1
        val = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        d *= val
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
        if PLOT1:
            if False:
                p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(d))
                p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(n.fft.fft(__d) * n.fft.fft(ker)))
                p.semilogy(freqs[MIN_CH:MAX_CH], n.abs(d - n.fft.fft(__d) * n.fft.fft(ker)))
                p.ylim(1e2,1e6)
            else:
                p.plot(freqs[MIN_CH:MAX_CH], d, '.')
                p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(__d) * n.fft.fft(ker), ',')
                #p.plot(freqs[MIN_CH:MAX_CH], d - n.fft.fft(__d) * n.fft.fft(ker))
        _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(ker))
        if True:
            ker = n.fft.ifft(val)
            __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            if DECONV1WIDTH == 0: __d[1:] = 0
            else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
            if PLOT1:
                #p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(__d) * n.fft.fft(ker))
                p.plot(freqs[MIN_CH:MAX_CH], n.fft.fft(_d) - n.fft.fft(__d) * n.fft.fft(ker), '.')
            _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(ker))
        d = n.fft.fft(_d) * val
        #if n.random.uniform() < .5: d = -d
        dat[i] = dat.get(i,0) + d
        wgt[i] = wgt.get(i,0) + val
    if PLOT1: p.show()

scores = {}
for i in dat:
    avg = dat[i] / wgt[i].clip(1,n.Inf)
    scores[i] = n.std(avg[100:-100])
mscore = n.median(scores.values())

for i in scores:
    print i, scores[i],
    if scores[i] > mscore:
        print '*'
        del(dat[i])
        del(wgt[i])
        continue
    print ''

if PLOT2:
    totavg = sum([dat[i] for i in dat]) / sum([wgt[i] for i in wgt]).clip(1,n.Inf)
    p.subplot(121)
    p.plot(freqs[MIN_CH:MAX_CH], avg, label='avg')
    p.subplot(122)
    p.plot(n.abs(n.fft.ifft(totavg[100:-100])), label='avg')
    p.legend()
    p.show()

data,gain = {},{}
__dat,__wgt,__dec,__cln = {},{},{},{}
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    for i in dat:
        _d = pfb(dat[i][ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
        _w = pfb(wgt[i][ch0:ch1], taps=NTAPS, window='hanning', fft=n.fft.ifft)
        g1 = n.abs(_w[0])
        if True:
            __d, info = a.deconv.clean(_d, _w, tol=opts.clean)
            __d[1:] = 0
            #__d[2:6] = 0; __d[7:] = 0
            _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(_w))
        _d1 = _d[:_d.size/2]
        for _d2,g2 in zip(data.get(b,[]), gain.get(b,[])):
            __dat[b] = __dat.get(b,0) + _d1 * n.conj(_d2)
            __wgt[b] = __wgt.get(b,0) + g1 * g2
        data[b] = data.get(b,[]) + [_d1]
        gain[b] = gain.get(b,[]) + [g1]

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
    k_pl = (eta * 1e-3 * dk_deta(z))[:__dat[b].size]
    #scalar = X2Y * k_pl**3 / (2*n.pi**2) * FoV * B / 8
    scalar = 1
    dat.append((__dat[b] / __wgt[b]) * scalar)
dat = n.array(dat).real
#avg = n.average(dat,axis=0); avg.shape = (1,avg.size); dat -= avg

#p.subplot(221)
dat1 = n.log10(n.abs(dat))
#mx,drng = 10,7
#mx,drng = 7,5
mx,drng = 3,4
#mx,drng = dat1.max(), dat1.max()-max(dat1.min(),2)
p.imshow(dat1, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
p.colorbar(shrink=.5)
print len(times)
p.show()


