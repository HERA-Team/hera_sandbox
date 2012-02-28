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
pb_poly = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]

def jy2T(freqs_in_GHz):
    lam = a.const.c / (freqs_in_GHz * 1e9)
    pb = n.polyval(pb_poly, freqs_in_GHz)
    return 1e-23 * lam**2 / (2 * a.const.k * pb)

def dk_deta(z):
    '''Mpc/MHz'''
    return 5.1 * 0.073 * ((1+z) / 10.)**-.5

def dL_dtheta(z):
    '''Mpc/radian'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

def dk_du(z):
    return 1. / (dL_dtheta(z) * (n.pi / 2))

def dL_dFQ(z):
    '''Mpc/MHz'''
    return 1.7 / 0.1 * ((1+z) / 10.)**.5

def early_eor_pspec(k):
    return k**-2 / 1e6

def late_eor_pspec(k):
    return 2*k**-3 / 1e6

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
            #freqs = f['freqs'][100:200]
            freqs = f['freqs'][50:150]
            metadata[fq]['freqs'] = freqs
            #metadata[fq]['dlys'] = n.arange(0, .5/(freqs[1]-freqs[0]) + 1./(freqs[-1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
            metadata[fq]['dlys'] = n.arange(0, 1./(freqs[1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
            metadata[fq]['dlys'] = n.where(metadata[fq]['dlys'] > .5/(freqs[1]-freqs[0]), metadata[fq]['dlys'] - (1./(freqs[1]-freqs[0])), metadata[fq]['dlys'])
            metadata[fq]['jy2T'] = n.average(jy2T(freqs))
            z = metadata[fq]['z'] = n.average(1.420405 / freqs - 1)
            L = metadata[fq]['L'] = freqs * 1e3 * dL_dFQ(z)
            #metadata[fq]['k_pl'] = n.arange(0, .5/(L[1]-L[0]) + 1.5/(L[-1]-L[0]), 1./(L[-1]-L[0]))
            metadata[fq]['k_pl'] = (metadata[fq]['dlys'] * 1e-3 * dk_deta(z))[:metadata[fq]['dlys'].size/2+1]
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
def fit_pwrlaw_noise(k, pspec, klow=.1, kns=1.):
    valid = n.where(k > klow, 1, 0)
    valid_noise = n.where(k > kns, 1, 0)
    k2 = k**2
    noise = n.sum((k2*pspec).compress(valid_noise)) / n.sum(k2.compress(valid_noise))
    noise_err = n.sqrt(n.sum((k2*pspec).compress(valid_noise)**2) / n.sum(k2.compress(valid_noise)))
    return noise, noise_err

#UMAG_RES = .2
#UMAG_RES = .25
UMAG_RES = 1.

print 'Sorting weights'
maxwgt = {}
for fq in _wgt:
    maxwgt[fq] = {}
    bins = _wgt[fq].keys()
    def cmp(a,b):
        if n.sum(n.abs(_wgt[fq][a])**2) > n.sum(n.abs(_wgt[fq][b])**2): return 1
        else: return -1
    bins.sort(cmp)
    for bin in bins:
        u,v,lst = bin2uv(bin,uv_res=UV_RES)
        #umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        umag = 2**(n.around(n.log2(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        maxwgt[fq][umag] = n.average(n.abs(_wgt[fq][bin])**2)

print 'Calculating Power Spectrum'
pspec_sum,pspec_wgt = {}, {}
pspec_mdl = {}
for fq in _sum:
    pspec_sum[fq],pspec_wgt[fq] = {},{}
    pspec_mdl[fq] = {}

for fq in _sum:
    for bin in _sum[fq]:
        wgt = n.average(n.abs(_wgt[fq][bin])**2)
        if wgt == 0: continue
        u,v,lst = bin2uv(bin,uv_res=UV_RES)
        #umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        umag = 2**(n.around(n.log2(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        #if wgt == 0 or wgt < WGT_FRAC * maxwgt[fq][umag]: continue
        avg_z = n.average(metadata[fq]['z'])
        k_pl = metadata[fq]['k_pl']
        k_pr = dk_du(avg_z) * umag
        ks = n.sqrt(k_pl**2 + k_pr**2)
        # Weights are already SNR**2, so no need to square again
        _d = metadata[fq]['jy2T'] * _sum[fq][bin]
        ker = _wgt[fq][bin]
        val = n.fft.fft(ker)
        gain = n.abs(ker[0])**2
        if True:
            _d,info = a.deconv.clean(_d, ker, tol=opts.clean)
            _d += info['res'] / gain
        else: _d /= gain
        # Generate a model of the smooth-spectrum sky for this data by extracting lower delay modes
        #sky_maxdly = umag / metadata[fq]['freqs'][-1]
        #if False:
        #    skymask = n.where(n.abs(metadata[fq]['dlys']) <  1.5*sky_maxdly, 1., 0)
        #    mdl = n.fft.fft(_d * skymask * val)
        #    _mdl = n.fft.ifft(mdl * win)
        #else:
        #    freqs = metadata[fq]['freqs']
        #    def tsky_sync(amp,index):
        #        return amp * (freqs / .150)**index
        #    def tsky_psrc(a, skymask, avg_index=-1.):
        #        amp = n.random.normal(size=skymask.size).astype(n.complex) * a * skymask
        #        phs = n.random.uniform(0,2*n.pi, size=skymask.size).astype(n.complex)
        #        return n.fft.ifft(amp * n.exp(1j*phs)) * n.sqrt(amp.size) * (freqs/.150)**avg_index
        #    NOISE_LEV = .00015
        #    #NOISE_LEV = .00030
        #    TPSRC_LEV = .6 * 1.5
        #    TSYNC_LEV = 1./200 * 1.5
        #    skymask = (1 - n.abs(metadata[fq]['dlys']) / (1.*sky_maxdly)).clip(0,1)
        #    Tsync = tsky_sync(240., -2.5)
        #    Tpsrc = tsky_psrc(TPSRC_LEV, skymask, -1)
        #    Tsky_none = Tsync * TSYNC_LEV + Tpsrc ; Tsky_none -= .9*n.average(Tsky_none)
        #    noise = NOISE_LEV * n.random.normal(scale=Tsync)
        #    mdl = Tsky_none + noise
        #    mdl *= val
        #    _mdl = n.fft.ifft(mdl)
        #if True:
        #    _mdl,info = a.deconv.clean(_mdl, ker, tol=opts.clean)
        #    _mdl += info['res'] / gain
        #else: _mdl /= gain
        _d *= n.sqrt(val.size)
        #_mdl *= n.sqrt(val.size)
        _pspec = n.abs(_d)**2
        _pspec[1:_pspec.size/2] += _pspec[:_pspec.size/2:-1]
        _pspec[1:_pspec.size/2] /= 2
        _pspec = _pspec[:_pspec.size/2+1]
        noise,noise_err = fit_pwrlaw_noise(n.abs(ks),_pspec)
        wgt = 1./noise_err**2
        # Weight each k sample by SNR**2 of power spectrum
        #pspec_sum[fq][umag] = pspec_sum[fq].get(umag,0) + n.abs(_d)**2 * wgt
        pspec_sum[fq][umag] = pspec_sum[fq].get(umag,0) + (n.abs(_d)**2 - noise) * wgt
        #pspec_mdl[fq][umag] = pspec_mdl[fq].get(umag,0) + n.abs(_mdl)**2 * wgt
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
    #m = pspec_mdl[fq][umag] / pspec_wgt[fq][umag]
    #if FOLD:
    #    m[1:m.size/2] += m[:m.size/2:-1]
    #    m[1:m.size/2] /= 2
    #    m = m[:m.size/2+1]
    #pspec_m[fq][umag] = m

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
            #(fa,fb,fc) = fit_pwrlaw_noise(_ks,_pspec)
            #p.loglog(_ks, n.abs(_pspec - fc), color+symbol1)
            p.loglog(_ks, _pspec, color+symbol1, label='%3.1f'%(avg_z))
            #_ks,_pspec_m = rebin_log(ks[3:], pspec_m[fq][umag][3:])
            #_ks = n.concatenate([ks[:3], _ks])
            #_pspec_m = n.concatenate([pspec_m[fq][umag][:3], _pspec_m])
            #(fa,fb,fc) = fit_pwrlaw_noise(_ks,_pspec_m)
            #p.loglog(_ks, n.abs(_pspec_m - fc), color+':')
            #p.loglog(_ks, _pspec_m, color+'.', label='%3.1f'%(avg_z))
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
    p.xlim(1e-2,2e0)
    #p.ylim(1e-3,1e0)
    p.ylim(1e-6,1e0)
    p.ylabel('$K^2$')
    p.xlabel(r'$k (h\ {\rm Mpc})^{-1}$')
    p.grid()

p.show()
