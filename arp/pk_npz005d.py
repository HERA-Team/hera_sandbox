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
o.add_option('--lst', dest='lsts', default='0:00.00_23:59.99')
o.add_option('--clean', dest='clean', type='float', default=1e-4)
opts, args = o.parse_args(sys.argv[1:])

lst_rng = map(lambda x: float(ephem.hours(x)), opts.lsts.split('_'))
pb_poly = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]

def uv2bin(u,v,lst,uv_res=1.5, lst_res=.2):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=1.5, lst_res=.2):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(LST_RES)
    return u,v, lst

def jy2T(freqs_in_GHz):
    lam = a.const.c / (freqs_in_GHz * 1e9)
    pb = n.polyval(pb_poly, freqs_in_GHz)
    return 1e-23 * lam**2 / (2 * a.const.k * pb)

def dk_deta(z):
    '''(1/Mpc)/(1/MHz)'''
    return 5.1 * 0.073 * ((1+z) / 10.)**-.5

def dL_dtheta(z):
    '''Mpc/radian'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

def dk_du(z):
    return 1. / (dL_dtheta(z) * (n.pi / 2))

def dL_dFQ(z):
    '''Mpc/MHz'''
    return 1.7 / 0.1 * ((1+z) / 10.)**.5

def rebin_log(x, y, nbins=10):
    '''For y=f(x), bin x into logrithmic (base 10) bins, and average y over
    these bin sizes.'''
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

SUBBAND = 50
UMAG_RES = 1.
PLOT_KCUBE = False
bins = {}
metadata = {}
pspec_sum,pspec_wgt = {}, {}
pspec_mdl = {}

for filename in args:
    print 'Loading', filename
    f = n.load(filename)
    fq = n.average(f['freqs'])
    if not bins.has_key(fq):
        bins[fq], metadata[fq] = {}, {}
        pspec_sum[fq],pspec_wgt[fq] = {},{}
        pspec_mdl[fq] = {}
    UV_RES = f['uvres']
    LST_RES = f['lstres']
    freqs = f['freqs'][SUBBAND/2:3*SUBBAND/2]
    metadata[fq]['freqs'] = freqs
    eta = n.arange(0, 1./(freqs[1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
    eta = n.where(eta > .5/(freqs[1]-freqs[0]), eta - (1./(freqs[1]-freqs[0])), eta)
    #metadata[fq]['dlys'] = n.arange(0, 1./(freqs[1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
    #metadata[fq]['dlys'] = n.where(metadata[fq]['dlys'] > .5/(freqs[1]-freqs[0]), metadata[fq]['dlys'] - (1./(freqs[1]-freqs[0])), metadata[fq]['dlys'])
    metadata[fq]['dlys'] = eta
    metadata[fq]['jy2T'] = n.average(jy2T(freqs))
    z = metadata[fq]['z'] = n.average(1.420405 / freqs - 1)
    L = metadata[fq]['L'] = freqs * 1e3 * dL_dFQ(z)
    metadata[fq]['k_pl'] = (eta * 1e-3 * dk_deta(z))[:eta.size/2+1]
    for kd in [_k for _k in f.files if _k.startswith('sum')]:
        kw = 'wgt'+kd[len('sum'):]
        bin = int(kd.split('_')[-1])
        u,v,lst = bin2uv(bin,uv_res=UV_RES,lst_res=LST_RES)
        # Only process LST bins specified on command line
        if lst_rng[0] < lst_rng[1] and (lst < lst_rng[0] or lst > lst_rng[1]): continue
        if lst_rng[1] < lst_rng[0] and (lst > lst_rng[0] and lst < lst_rng[1]): continue
        _d1,_w1 = f[kd], f[kw]
        _d1 *= metadata[fq]['jy2T']   # Put into temp units
        # When we divide by g1 below, it contains a factor of size, and has the same FFT normalization as _d1.
        # Therefore, dividing by g1 puts the result back into temp units, but divides temp units by a factor of size.
        # Parseval's Thm says that temp^2 should be divided by factor of size, so we over-divided by sqrt(size).
        # This next step fixes that.
        #_d1 *= n.sqrt(_w1.size) # wrong
        # Parseval's Thm says sig_i^2 = sum(sig_k^2).  n.fft.ifft does that.
        # Dividing by g1 then removes any additional weighting
        # So upshot is no normalization beyond /g1 is required
        g1 = n.abs(_w1[0])
        if g1 == 0: continue
        if True:    # Deconvolve point sources
            _d1,info = a.deconv.clean(_d1, _w1, tol=opts.clean)
            _d1 += info['res'] / g1
        else: _d1 /= g1
        # Cross-corr this sample with all others for this bin, avoiding
        # Auto-corrs that introduce a noise^2 into P(k)
        for (_d2,g2) in bins[fq].get(bin,[]):
            g12 = g1 * g2
            pspec_sum[fq][bin] = pspec_sum[fq].get(bin,0) + _d1*n.conj(_d2) * g12
            pspec_wgt[fq][bin] = pspec_wgt[fq].get(bin,0) + g12
        bins[fq][bin] = bins[fq].get(bin,[]) + [(_d1,g1)]
        #print ephem.hours(lst), len(bins[fq][bin])

pkumag_sum,pkumag_wgt = {},{}
for fq in pspec_sum:
    pkumag_sum[fq],pkumag_wgt[fq] = {},{}
    for bin in pspec_sum[fq]:
        _d2,_w2 = pspec_sum[fq][bin], pspec_wgt[fq][bin]
        u,v,lst = bin2uv(bin,uv_res=UV_RES)
        umag = 2**(n.around(n.log2(n.sqrt(u**2+v**2)) / UMAG_RES) * UMAG_RES)
        avg_z = n.average(metadata[fq]['z'])
        k_pl = metadata[fq]['k_pl']
        k_pr = dk_du(avg_z) * umag
        ks = n.sqrt(k_pl**2 + k_pr**2)
        # Weight each k sample by SNR**2 of power spectrum
        pkumag_sum[fq][umag] = pkumag_sum[fq].get(umag,0) + _d2
        pkumag_wgt[fq][umag] = pkumag_wgt[fq].get(umag,0) + _w2

pspec = {}
pspec_m = {}
for fq in pspec_sum:
  pspec[fq] = {}
  pspec_m[fq] = {}
  for umag in pkumag_sum[fq]:
    d = n.abs(pkumag_sum[fq][umag] / pkumag_wgt[fq][umag])
    d[1:d.size/2] += d[:d.size/2:-1]
    d[1:d.size/2] /= 2
    d = d[:d.size/2+1]
    pspec[fq][umag] = d

nplots = {}
for fq in pspec:
    for umag in pspec[fq]:
        nplots[umag] = None
umags = nplots.keys(); umags.sort()
nplots = len(umags)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)

fqs = pspec.keys(); fqs.sort()

print 'Generating output files'
for cnt, umag in enumerate(umags):
    p.subplot(d1,d2,cnt + 1)
    for i,fq in enumerate(fqs):
        if not pspec[fq].has_key(umag): continue
        avg_z = n.average(metadata[fq]['z'])
        k_pl = metadata[fq]['k_pl']
        k_pr = dk_du(avg_z) * umag
        ks = n.sqrt(k_pl**2 + k_pr**2)
        _ks,_pspec = rebin_log(ks[3:], pspec[fq][umag][3:])
        _ks = n.concatenate([ks[:3], _ks])
        _pspec = n.concatenate([pspec[fq][umag][:3], _pspec])
        filename = 'pk_npz__%03d__%5.3f.npz' % (umag,fq)
        print 'Writing', filename
        n.savez(filename, k_pl=k_pl, k_pr=k_pr, ks=ks, _ks=_ks, pspec=pspec[fq][umag], _pspec=_pspec)
