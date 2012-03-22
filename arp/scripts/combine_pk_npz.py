#!/usr/bin/env python
import numpy as n
import os, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('pk_npz.py [options] *.npz')
opts, args = o.parse_args(sys.argv[1:])

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
            metadata[fq]['jy2T'] = jy2T(freqs)
            z = metadata[fq]['z'] = 1.420405 / freqs - 1
            L = metadata[fq]['L'] = freqs * 1e3 * dL_dFQ(z)
            metadata[fq]['k_pl'] = n.arange(0, .5/(L[1]-L[0]) + 1./(L[-1]-L[0]), 1./(L[-1]-L[0]))
        elif k.startswith('sum'):
            bin = int(k.split('_')[-1])
            _sum[fq][bin] = _sum[fq].get(bin,0) + f[k]
        elif k.startswith('wgt'):
            bin = int(k.split('_')[-1])
            _wgt[fq][bin] = _wgt[fq].get(bin,0) + f[k]
        else: raise ValueError('Unrecognized keyword: %s' % k)

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
        u,v = bin2uv(bin,uv_res=UV_RES)
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / .2) * .2)
        maxwgt[fq][umag] = n.average(_wgt[fq][bin]**2)

print 'Calculating Power Spectrum'
pspec_sum,pspec_wgt = {}, {}
for fq in _sum: pspec_sum[fq],pspec_wgt[fq] = {},{}

for fq in _sum:
    for bin in _sum[fq]:
        wgt = n.average(_wgt[fq][bin]**2)
        u,v = bin2uv(bin,uv_res=UV_RES)
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / .2) * .2)
        if wgt < WGT_FRAC * maxwgt[fq][umag]: continue
        d = metadata[fq]['jy2T'] * _sum[fq][bin] / n.clip(_wgt[fq][bin], 1, n.Inf)
        val = n.where(_wgt[fq][bin] < _wgt[fq][bin].max()/2, 0, 1.)
        d *= val
        # Add some windowing to reduce sidelobes of smooth-spec sky
        if False: win = window_hanning(val.size)
        elif False: win = window_blackman_harris(val.size)
        else: win = n.ones(val.size)
        _d = n.fft.ifft(d*win)
        ker = n.fft.ifft(val*win) 
        gain = n.abs(ker[0])**2
        if not n.any(val): continue
        if True:
            _d,info = a.deconv.clean(_d, ker, tol=1e-4)
            _d += info['res'] / gain
        else: _d /= gain
        d = n.fft.fft(_d)
        _d *= n.sqrt(val.size)
        # Weight each k sample by SNR**2 of power spectrum
        pspec_sum[fq][umag] = pspec_sum[fq].get(umag,0) + n.abs(_d)**2 * wgt
        pspec_wgt[fq][umag] = pspec_wgt[fq].get(umag,0) + wgt

pspec = {}
for fq in pspec_sum:
  pspec[fq] = {}
  for umag in pspec_sum[fq]:
    d = pspec_sum[fq][umag] / pspec_wgt[fq][umag]
    d[1:d.size/2] += d[:d.size/2:-1]
    d[1:d.size/2] /= 2
    d = d[:d.size/2+1]
    pspec[fq][umag] = n.sqrt(d)

nplots = {}
for fq in pspec:
    for umag in pspec[fq]:
        nplots[umag] = None
umags = nplots.keys(); umags.sort()
nplots = len(umags)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)
#p.subplot(d1,d2,1)
#for i,fq in enumerate(fqs):
#    color = 'kbgrcmy'[i%7]
#    if i >= 7: symbol = '--'
#    else: symbol = '-'
#    for umag in pspec[fq]:
#        z = n.average(metadata[fq]['z'])
#        p.semilogy(pspec[fq][umag], color+symbol, label='z=%3.1f%s'%(z,k))
#        p.title(str(umag))
#        break
##p.legend()
#p.ylim(1e-3,1e0)
#p.ylabel('K')
#p.xlabel('$k_\|$')
#p.grid()

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
        if i >= 7: symbol = '--'
        else: symbol = '-'
        _ks,_pspec = rebin_log(ks[3:], pspec[fq][umag][3:])
        _ks = n.concatenate([ks[:3], _ks])
        _pspec = n.concatenate([pspec[fq][umag][:3], _pspec])
        p.loglog(_ks, _pspec, color+symbol, label='%3.1f'%(avg_z))
    p.title(str(umag))
    p.xlim(1e-2,1e1)
    p.ylim(1e-3,1e0)
    p.ylabel('K')
    p.xlabel('$k (h\ {\rm Mpc})^{-1}$')
    p.grid()

p.show()
