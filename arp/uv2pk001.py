#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('filter_src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
src = a.phs.RadioFixedBody(ra='12:30', dec='40:00')

inttime = uv['inttime']
MIN_CH,MAX_CH,SUBBAND = 150,690,60
#MIN_CH,MAX_CH,SUBBAND = 150,690,36
#MIN_CH,MAX_CH,SUBBAND = 240,720,120
aa.select_chans(n.arange(MIN_CH, MAX_CH))
#UV_RES = 1.5
UV_RES = 0.5
WGT_FRAC = .5
pb_poly_old = [1.02854332e+09, -9.49707493e+08, 3.64775002e+08, -7.46038156e+07, 8.56951433e+06, -5.24246222e+05, 1.33464786e+04]
pb_poly = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06,
    7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]

def jy2T(freqs_in_GHz):
    lam = a.const.c / (freqs_in_GHz * 1e9)
    pb = n.polyval(pb_poly, freqs_in_GHz)
    return 1e-23 * lam**2 / (2 * a.const.k * pb)

def uv2bin(u,v, uv_res=UV_RES):
    return int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)
def bin2uv(bin, uv_res=UV_RES):
    v = (bin % 8192 - 4096) * float(uv_res)
    u = (bin / 8192 - 4096) * float(uv_res)
    return u,v

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

freqs = aa.get_afreqs()
jy2T = jy2T(freqs)
zs = 1420.05 / freqs - 1
Ls = freqs * dL_dFQ(zs)
ks = n.arange(0, .5/(Ls[1]-Ls[0]) + 1./(Ls[-1]-Ls[0]), 1./(Ls[-1]-Ls[0]))
ks = ks.clip(ks[1]/4, n.Inf)


_sum,_wgt = {},{}
for k in 'abc': _sum[k],_wgt[k] = {}, {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    times = []
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            aa.set_jultime(t)
            src.compute(aa)
            times.append(t)
        #if len(times) < 20 or len(times) > 230: continue
        bl = a.miriad.ij2bl(i,j)
        u,v,w = aa.get_baseline(i,j, src=src)
        u *= aa.get_afreqs().max()
        v *= aa.get_afreqs().max()
        if n.abs(u) < 12.5: continue # overly short baselines adversely affected by xtalk removal
        d = d[MIN_CH:MAX_CH]
        f = f[MIN_CH:MAX_CH]
        d = aa.phs2src(d, src, i, j)
        d /= aa.passband(i,j)
        #print n.median(n.abs(d)),
        val = n.logical_not(f).astype(n.float)
        d *= jy2T * val
        #print n.median(n.abs(d))
        if not n.all(val == 0):
            # Maybe do this more than once?
            if False:
                gain = n.sqrt(n.average(val**2))
                ker = n.fft.ifft(val)
                _d = n.fft.ifft(d)
                _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
                if True: _d = info['res'] / gain
                else: _d += info['res'] / gain
                d = n.fft.fft(_d) * val
            # accumulate for different amounts of time
            for cnt,k in enumerate('abc'):
                if not _sum[k].has_key(bl): _sum[k][bl],_wgt[k][bl] = {},{}
                # ToDo: consider that lower freqs should change bins less often
                # than higher freqs, make u,v freq-dependent
                bin = uv2bin(u,v,uv_res=UV_RES*(cnt+1))
                #_sum[k][bl][bin] = _sum[k][bl].get(bin,0) + n.fft.fft(_d) * val / n.clip(aa.passband(i,j)*val, 4e-6, n.Inf) * jy2T
                _sum[k][bl][bin] = _sum[k][bl].get(bin,0) + d
                _wgt[k][bl][bin] = _wgt[k][bl].get(bin,0) + val

#for k in _sum:
#    for bl in _sum[k]:
#        for bin in _sum[k][bl]:
#            d = _sum[k][bl][bin] / _wgt[k][bl][bin].clip(1,n.Inf)
#            print k, a.miriad.bl2ij(bl), bin, n.median(n.abs(d[1:] - d[:-1])),
#            print n.median(n.abs(d[1:] - d[:-1]) * n.sqrt(_wgt[k][bl][bin][1:]))

maxwgt = {}
for k in _wgt:
    maxwgt[k] = {}
    blbins = [(bl,bin) for bl in _wgt[k] for bin in _wgt[k][bl]]
    def cmp(a,b):
        abl,abin = a
        bbl,bbin = b
        if n.sum(_wgt[k][abl][abin]) > n.sum(_wgt[k][bbl][bbin]): return 1
        else: return -1
    blbins.sort(cmp)
    for (bl,bin) in blbins:
        u,v = bin2uv(bin,uv_res=UV_RES*('abcde'.find(k)+1))
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / .2) * .2)
        maxwgt[k][umag] = n.average(_wgt[k][bl][bin]**2)
    #bl,bin = blbins[-1]
    #for (bl,bin) in blbins:
    #    if n.average(_wgt[k][bl][bin]**2) < WGT_FRAC * maxwgt[k]: continue
    #    print k, a.miriad.bl2ij(bl), bin2uv(bin,uv_res=UV_RES*('abcde'.find(k)+1)),
    #    print n.around(n.sqrt(n.average(_wgt[k][bl][bin]**2))*inttime)

import pylab as p
pspec_sum,pspec_wgt = {}, {}
ker_sum = {}
for k in _sum:
    pspec_sum[k],pspec_wgt[k] = {},{}
    ker_sum[k] = {}

for k in _sum:
  for bl in _sum[k]:
    i,j = a.miriad.bl2ij(bl)
    for bin in _sum[k][bl]:
        wgt = n.average(_wgt[k][bl][bin]**2)
        u,v = bin2uv(bin,uv_res=UV_RES*('abcde'.find(k)+1))
        umag = 10**(n.around(n.log10(n.sqrt(u**2+v**2)) / .2) * .2)
        if wgt < WGT_FRAC * maxwgt[k][umag]: continue
        d = _sum[k][bl][bin] / n.clip(_wgt[k][bl][bin], 1, n.Inf)
        val = n.where(_wgt[k][bl][bin] < _wgt[k][bl][bin].max()/2, 0, 1.)
        d *= val
        d.shape = (d.size/SUBBAND,SUBBAND)
        _wgt[k][bl][bin].shape = (d.size/SUBBAND,SUBBAND)
        val.shape = (val.size/SUBBAND,SUBBAND)
        # Add some windowing to reduce sidelobes of smooth-spec sky
        if True: win = window_hanning(val.shape[1])
        elif False: win = window_blackman_harris(val.shape[1])
        else: win = n.ones(val.shape[1])
        win.shape = (1,win.size)
        _d = n.fft.ifft(d*win, axis=1)
        #gain = n.sqrt(n.average(val**2, axis=1))
        ker = n.fft.ifft(val*win, axis=1) 
        gain = n.abs(ker[:,0])**2
        for cnt in range(d.shape[0]):
            if not n.any(val[cnt]): continue
            if True:
                _d[cnt],info = a.deconv.clean(_d[cnt], ker[cnt], tol=1e-4)
                _d[cnt] += info['res'] / gain[cnt]
            else: _d[cnt] /= gain[cnt]
            d = n.fft.fft(_d[cnt])
            _d[cnt] *= n.sqrt(val[cnt].size)
            #print k, a.miriad.bl2ij(bl), bin, cnt, n.median(n.abs(d[1:] - d[:-1])),
            #print n.median(n.abs(d[1:] - d[:-1]) * n.sqrt(_wgt[k][bl][bin][cnt,1:])),
            #print n.median(n.abs(_d[cnt][-10:]))
        # Weight each k sample by SNR**2 of power spectrum
        pspec_sum[k][umag] = pspec_sum[k].get(umag,0) + n.abs(_d)**2 * wgt
        pspec_wgt[k][umag] = pspec_wgt[k].get(umag,0) + wgt
        #ker_sum[k][umag] = ker_sum[k].get(umag,0) + n.abs(ker)**2 * wgt

#for k in pspec_sum:
#    for umag in pspec_sum[k]:
#        d = n.pspec_sum[k][umag] / pspec_wgt[k][umag].clip(1,n.Inf)
#        print k, a.miriad.bl2ij(bl), bin, n.median(n.abs(d[1:] - d[:-1])),
#        print n.median(n.abs(d[1:] - d[:-1]) * n.sqrt(_wgt[k][bl][bin][1:]))

#for umag in ker_sum['c']:
#  print umag
#  for i in range(ker_sum[k][umag].shape[0]):
#    color = 'kbgrcmy'[i%7]
#    if i>= 7: symbol = '--'
#    else: symbol = '-'
#    ker_avg = ker_sum['c'][umag][i,1:ker_sum['c'][umag].shape[1]/2]
#    ker_avg = n.sqrt(ker_avg / pspec_wgt['c'][umag])
#    p.plot(ker_avg, color+symbol)
#  p.show()

pspec = {}
for k in pspec_sum:
  pspec[k] = {}
  for umag in pspec_sum[k]:
    #pspec_sum[k][umag] = pspec_sum[k][umag][:,1:pspec_sum[k][umag].shape[1]/2] + n.fliplr(pspec_sum[k][umag][:,pspec_sum[k][umag].shape[1]/2+1:])
    #pspec_wgt[k][umag] *= 2
    d = pspec_sum[k][umag] / pspec_wgt[k][umag]
    d[:,1:d.shape[1]/2] += n.fliplr(d[:,d.shape[1]/2+1:])
    d[:,1:d.shape[1]/2] /= 2
    d = d[:,:d.shape[1]/2+1]
    pspec[k][umag] = n.sqrt(d)


nplots = len(pspec.values()[0]) + 1
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)
p.subplot(d1,d2,1)
for k in pspec:
 if k == 'b': continue
 for umag in pspec[k]:
  for i in range(pspec[k][umag].shape[0]):
    fq = freqs[SUBBAND*i+SUBBAND/2]
    z = n.around(1.420405 / fq - 1, 1)
    #print SUBBAND*i+SUBBAND/2, fq, z
    color = 'kbgrcmy'[i%7]
    if i >= 7: symbol = '--'
    else: symbol = '-'
    p.semilogy(pspec[k][umag][i], 
        color+symbol, label='z=%3.1f%s'%(z,k))
  p.title(str(umag))
  break
#p.legend()
p.ylim(1e-3,1e0)
p.ylabel('K')
p.xlabel('$k_\|$')
p.grid()

def rebin_log(x, y, nbins=10):
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, range=(logx[1],logx[-1]), bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, range=(logx[1],logx[-1]), bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

for k in pspec:
  for cnt, umag in enumerate(pspec[k]):
    if k != 'c': continue
    #p.subplot(222 + 'abcde'.find(k))
    p.subplot(d1,d2,cnt + 2)
    for i in range(pspec[k][umag].shape[0]):
        fq = freqs[SUBBAND*i:SUBBAND*(i+1)]
        zs = 1.420405 / fq - 1
        Ls = 1e3*fq * dL_dFQ(zs)
        k_pl = n.arange(0, .5/(Ls[1]-Ls[0]) + 1./(Ls[-1]-Ls[0]), 1./(Ls[-1]-Ls[0]))
        k_pr = dk_du(n.average(zs)) * umag
        ks = n.sqrt(k_pl**2 + k_pr**2)
        color = 'kbgrcmy'[i%7]
        if i >= 7: symbol = '--'
        else: symbol = '-'
        _ks,_pspec = rebin_log(ks, pspec[k][umag][i])
        #p.loglog(ks, pspec[k][umag][i], color+symbol, label='%3.1f'%(n.average(zs)))
        p.loglog(_ks, _pspec, color+symbol, label='%3.1f'%(n.average(zs)))
    if k == 'b': p.legend()
    p.title(str(umag))
    p.xlim(1e-2,1e1)
    p.ylim(1e-3,1e0)
    p.ylabel('K')
    p.xlabel('$k_\|$')
    p.grid()

p.show()
