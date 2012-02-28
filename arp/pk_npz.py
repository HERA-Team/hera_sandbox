#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, pylab as p
import os, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('pk_npz.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--clean', dest='clean', type='float', default=1e-4)
opts, args = o.parse_args(sys.argv[1:])

WGT_FRAC = .5
#pb_poly_old = [1.02854332e+09, -9.49707493e+08, 3.64775002e+08, -7.46038156e+07, 8.56951433e+06, -5.24246222e+05, 1.33464786e+04]
pb_poly = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]

def jy2T(freqs_in_GHz):
    lam = a.const.c / (freqs_in_GHz * 1e9)
    pb = n.polyval(pb_poly, freqs_in_GHz)
    return 1e-23 * lam**2 / (2 * a.const.k * pb)

def window(N):
    return 0.35875 - 0.48829 * n.cos(2*n.pi*n.arange(N)/N) + 0.14128 * n.cos(4*n.pi*n.arange(N)/N) - 0.01168 * n.cos(6*n.pi *n.arange(N)/N)

def dL_dFQ(z):
    '''Mpc/MHz'''
    return 1.7 / 0.1 * ((1+z) / 10.)**.5

def dL_dtheta(z):
    '''Mpc/radian'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

def dk_du(z):
    return 1. / (dL_dtheta(z) * (n.pi / 2))

def early_eor_pspec(k):
    return k**-2 / 1e6

def late_eor_pspec(k):
    return 2*k**-3 / 1e6

def window_blackman_harris(N):
    '''Blackman-Harris -- minimizes sidelobes at major SNR expense'''
    x = n.arange(N,dtype=n.float) / (N-1)
    return 0.35875 - 0.48829 * n.cos(2*n.pi*x) + 0.14128 * n.cos(4*n.pi*x) - 0.01168 * n.cos(6*n.pi*x)

def window_hanning(N):
    x = n.arange(N,dtype=n.float) / (N-1)
    return 0.5 * (1 - n.cos(2*n.pi*x))

_sum,_wgt = {},{}
metadata = {}
for filename in args:
    print 'Loading', filename
    f = n.load(filename)
    fq = n.average(f['freqs'])
    if not _sum.has_key(fq):
        _sum[fq],_wgt[fq] = {},{}
        metadata[fq] = {}
    for k in f.files:
        if k == 'freqs':
            freqs = f['freqs']
            metadata[fq]['freqs'] = freqs
            #metadata[fq]['dlys'] = n.arange(0, .5/(freqs[1]-freqs[0]) + 1./(freqs[-1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
            metadata[fq]['dlys'] = n.arange(0, 1./(freqs[1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
            metadata[fq]['dlys'] = n.where(metadata[fq]['dlys'] > .5/(freqs[1]-freqs[0]), metadata[fq]['dlys'] - (1./(freqs[1]-freqs[0])), metadata[fq]['dlys'])
            metadata[fq]['jy2T'] = jy2T(freqs)
            z = metadata[fq]['z'] = 1.420405 / freqs - 1
            L = metadata[fq]['L'] = freqs * 1e3 * dL_dFQ(z)
            metadata[fq]['k_pl'] = n.arange(0, .5/(L[1]-L[0]) + 1.5/(L[-1]-L[0]), 1./(L[-1]-L[0]))
        elif k == 'uvres':
            UV_RES = f['uvres']
        elif k == 'lstres':
            LST_RES = f['lstres']
        elif k.startswith('sum'):
            bin = int(k.split('_')[-1])
            _sum[fq][bin] = _sum[fq].get(bin,0) + f[k]
        elif k.startswith('wgt'):
            bin = int(k.split('_')[-1])
            _wgt[fq][bin] = _wgt[fq].get(bin,0) + f[k]
        else: raise ValueError('Unrecognized keyword: %s' % k)

def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(LST_RES)
    return u,v, lst

#UMAG_RES = .2
UMAG_RES = .25

print 'Sorting weights'
maxwgt = {}
for fq in _wgt:
    maxwgt[fq] = {}
    bins = _wgt[fq].keys()
    def cmp(a,b):
        if n.sum(_wgt[fq][a]) > n.sum(_wgt[fq][b]): return 1
        else: return -1
    bins.sort(cmp)
    for bin in bins:
        u,v,lst = bin2uv(bin,uv_res=UV_RES)
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        maxwgt[fq][umag] = n.average(_wgt[fq][bin]**2)

print 'Calculating Power Spectrum'
pspec_sum,pspec_wgt = {}, {}
pspec_mdl = {}
for fq in _sum:
    pspec_sum[fq],pspec_wgt[fq] = {},{}
    pspec_mdl[fq] = {}

for fq in _sum:
    for bin in _sum[fq]:
        wgt = n.average(_wgt[fq][bin]**2)
        u,v,lst = bin2uv(bin,uv_res=UV_RES)
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        if wgt < WGT_FRAC * maxwgt[fq][umag] or wgt == 0: continue
        # Add some windowing to reduce sidelobes of smooth-spec sky
        #if True: win = window_hanning(_wgt[fq][bin].size)
        if True: win = window_blackman_harris(_wgt[fq][bin].size)
        else: win = n.ones(_wgt[fq][bin].size)
        if False:
            d = metadata[fq]['jy2T'] * _sum[fq][bin] / n.clip(_wgt[fq][bin], 1, n.Inf)
            val = n.where(_wgt[fq][bin] < _wgt[fq][bin].max()/2, 0, 1.)
            d *= val
            _d = n.fft.ifft(d*win)
            ker = n.fft.ifft(val*win) 
        else:
            d = metadata[fq]['jy2T'] * _sum[fq][bin]
            val = _wgt[fq][bin]
            # Weights are already SNR**2, so no need to square again
            _d = n.fft.ifft(d*win)
            ker = n.fft.ifft(val*win) 
        gain = n.abs(ker[0])**2
        if True:
            _d,info = a.deconv.clean(_d, ker, tol=opts.clean)
            _d += info['res'] / gain
        else: _d /= gain
        # Generate a model of the smooth-spectrum sky for this data by extracting lower delay modes
        sky_maxdly = umag / metadata[fq]['freqs'][-1]
        if False:
            skymask = n.where(n.abs(metadata[fq]['dlys']) <  1.5*sky_maxdly, 1., 0)
            mdl = n.fft.fft(_d * skymask * val)
            _mdl = n.fft.ifft(mdl * win)
        else:
            freqs = metadata[fq]['freqs']
            def tsky_sync(amp,index):
                return amp * (freqs / .150)**index
            def tsky_psrc(a, skymask, avg_index=-1.):
                amp = n.random.normal(size=skymask.size).astype(n.complex) * a * skymask
                phs = n.random.uniform(0,2*n.pi, size=skymask.size).astype(n.complex)
                return n.fft.ifft(amp * n.exp(1j*phs)) * n.sqrt(amp.size) * (freqs/.150)**avg_index
            NOISE_LEV = .00015
            #NOISE_LEV = .00030
            TPSRC_LEV = .6 * 1.5
            TSYNC_LEV = 1./200 * 1.5
            skymask = (1 - n.abs(metadata[fq]['dlys']) / (1.*sky_maxdly)).clip(0,1)
            Tsync = tsky_sync(240., -2.5)
            Tpsrc = tsky_psrc(TPSRC_LEV, skymask, -1)
            Tsky_none = Tsync * TSYNC_LEV + Tpsrc ; Tsky_none -= .9*n.average(Tsky_none)
            noise = NOISE_LEV * n.random.normal(scale=Tsync)
            mdl = Tsky_none + noise
            mdl *= val
            _mdl = n.fft.ifft(mdl * win)
        if True:
            _mdl,info = a.deconv.clean(_mdl, ker, tol=opts.clean)
            _mdl += info['res'] / gain
        else: _mdl /= gain
        _d *= n.sqrt(val.size)
        _mdl *= n.sqrt(val.size)
        # Weight each k sample by SNR**2 of power spectrum
        pspec_sum[fq][umag] = pspec_sum[fq].get(umag,0) + n.abs(_d)**2 * wgt
        pspec_mdl[fq][umag] = pspec_mdl[fq].get(umag,0) + n.abs(_mdl)**2 * wgt
        pspec_wgt[fq][umag] = pspec_wgt[fq].get(umag,0) + wgt

FOLD = True
pspec = {}
pspec_m = {}
for fq in pspec_sum:
  pspec[fq] = {}
  pspec_m[fq] = {}
  for umag in pspec_sum[fq]:
    d = pspec_sum[fq][umag] / pspec_wgt[fq][umag]
    if FOLD:
        d[1:d.size/2] += d[:d.size/2:-1]
        d[1:d.size/2] /= 2
        d = d[:d.size/2+1]
    pspec[fq][umag] = d
    m = pspec_mdl[fq][umag] / pspec_wgt[fq][umag]
    if FOLD:
        m[1:m.size/2] += m[:m.size/2:-1]
        m[1:m.size/2] /= 2
        m = m[:m.size/2+1]
    pspec_m[fq][umag] = m

nplots = {}
for fq in pspec:
    for umag in pspec[fq]:
        nplots[umag] = None
umags = nplots.keys(); umags.sort()
nplots = len(umags)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)

def rebin_log(x, y, nbins=10):
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

fqs = pspec.keys(); fqs.sort()

def fit_pwrlaw_noise(k, pspec, klow=.03, kns=.4):
    valid = n.where(k > klow, 1, 0)
    valid_noise = n.where(k > kns, 1, 0)
    k2 = k**2
    noise = n.sum((k2*pspec).compress(valid_noise)) / n.sum(k2.compress(valid_noise))
    #noise = n.min(pspec.compress(valid_noise))
    logk = n.log10(k.compress(valid))
    logp = n.log10(n.abs(pspec.compress(valid) - noise))
    poly = n.polyfit(logk, logp, deg=1)
    amp = 10**poly[-1]
    slope = poly[-2]
    prms = (amp,slope,noise)
    return prms

for cnt, umag in enumerate(umags):
    p.subplot(d1,d2,cnt + 1)
    for i,fq in enumerate(fqs):
        if not pspec[fq].has_key(umag): continue
        avg_z = n.average(metadata[fq]['z'])
        k_pl = metadata[fq]['k_pl']
        k_pr = dk_du(avg_z) * umag
        ks = n.sqrt(k_pl**2 + k_pr**2)
        color = 'kbgrcmy'[i%7]
        if i >= 7: symbol1,symbol2 = '--','^'
        else: symbol1,symbol2 = '-','v'
        if FOLD:
            _ks,_pspec = rebin_log(ks[3:], pspec[fq][umag][3:])
            _ks = n.concatenate([ks[:3], _ks])
            _pspec = n.concatenate([pspec[fq][umag][:3], _pspec])
            (fa,fb,fc) = fit_pwrlaw_noise(_ks,_pspec)
            #fa,fc = max(fa,1e-10), max(fc,1e-10)
            #p.loglog(_ks, fa*_ks**fb+fc, color+':')
            #p.loglog(_ks, fc*n.ones_like(_ks), color+':')
            p.loglog(_ks, n.abs(_pspec - fc), color+symbol2)
            #p.loglog(_ks, fa*_ks**fb, color+':')
            p.loglog(_ks, _pspec, color+symbol1, label='%3.1f'%(avg_z))
            _ks,_pspec_m = rebin_log(ks[3:], pspec_m[fq][umag][3:])
            _ks = n.concatenate([ks[:3], _ks])
            _pspec_m = n.concatenate([pspec_m[fq][umag][:3], _pspec_m])
            (fa,fb,fc) = fit_pwrlaw_noise(_ks,_pspec_m)
            p.loglog(_ks, n.abs(_pspec_m - fc), color+'.')
            p.loglog(_ks, _pspec_m, color+':', label='%3.1f'%(avg_z))
        else:
            # plot first half
            d = pspec[fq][umag]
            _ks,_pspec = rebin_log(ks[3:], d[3:d.size/2+1])
            _ks = n.concatenate([ks[:3], _ks])
            _pspec = n.concatenate([d[:3], _pspec])
            p.loglog(_ks, _pspec, color+symbol1, label='%3.1f'%(avg_z))
            # plot second half
            _ks,_pspec = rebin_log(ks[3:], d[:d.size/2-1:-1][2:])
            _ks = n.concatenate([ks[:3], _ks])
            _pspec = n.concatenate([d[0:1], d[:-3:-1], _pspec])
            p.loglog(_ks, _pspec, color+symbol1, label='%3.1f'%(avg_z))
    p.title(str(umag))
    #p.xlim(1e-2,1e1)
    p.xlim(1e-2,1e0)
    #p.ylim(1e-3,1e0)
    p.ylim(1e-6,1e0)
    p.ylabel('$K^2$')
    p.xlabel(r'$k (h\ {\rm Mpc})^{-1}$')
    p.grid()

p.show()
