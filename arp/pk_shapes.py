#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a

def tsky_sync(amp,index):
    return amp * (freqs / .150)**index

def dL_dFQ(z):
    '''Mpc/MHz'''
    return 1.7 / 0.1 * ((1+z) / 10.)**.5

def dL_dtheta(z):
    '''Mpc/radian'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2

def dk_du(z):
    #return 1. / (dL_dtheta(z) * (n.pi / 2))
    return 2*n.pi / (dL_dtheta(z) * (n.pi / 2))

def early_eor_pspec(k):
    return k**-2 / 1e6

def late_eor_pspec(k):
    return n.where(k >= .1, 2*k**-3 / 1e6, 2*(.1)**-3/1e6 * (k/.1)**-1)

def pspec2spec(pspec):
    amp = n.random.normal(scale=n.sqrt(pspec)).astype(n.complex)
    phs = n.random.uniform(0,2*n.pi, size=amp.size).astype(n.complex)
    return n.fft.irfft(amp * n.exp(1j*phs)) * n.sqrt(amp.size)

def window_blackman_harris(N):
    '''Blackman-Harris -- minimizes sidelobes at major SNR expense'''
    x = n.arange(N,dtype=n.float) / (N-1)
    return 0.35875 - 0.48829 * n.cos(2*n.pi*x) + 0.14128 * n.cos(4*n.pi*x) - 0.01168 * n.cos(6*n.pi*x)

def window_hanning(N):
    x = n.arange(N,dtype=n.float) / (N-1)
    return 0.5 * (1 - n.cos(2*n.pi*x))

def tsky_psrc(a, skymask, avg_index=-1.):
    amp = n.random.normal(size=skymask.size).astype(n.complex) * a * skymask
    phs = n.random.uniform(0,2*n.pi, size=skymask.size).astype(n.complex)
    return n.fft.irfft(amp * n.exp(1j*phs)) * n.sqrt(amp.size) * (freqs/.150)**avg_index

def rebin_log(x, y, nbins=10):
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

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

#UMAGS = [17.78, 31.62, 56.23, 100.00, 177.83]
#UMAGS = [25., 50., 100., 200.]
UMAGS = [16., 32., 64., 128.]
nplots = len(UMAGS)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)

for ucnt, UMAG in enumerate(UMAGS):
  p.subplot(d1,d2,ucnt + 1)
  #for cnt, FQ1 in enumerate([.134,.141,.148,.156,.159,.167,.172]):
  for cnt, FQ1 in enumerate([.134,.145,.156,.167]):
    #FQ1,FQ2,dFQ = .150, .15599, .0001
    FQ2,dFQ = FQ1+.00599, .0001
    #UMAG = [15.85, 25.12, 39.81, 63.10, 100.00, 158.49][2]
    #UMAG = [17.78, 31.62, 56.23, 100.00, 177.83][2]
    #skymax_kbin = UMAG / FQ2 * (FQ2 - FQ1)
    freqs = n.arange(FQ1,FQ2,dFQ)
    dlys = n.arange(0, .5/(freqs[1]-freqs[0]) + 1./(freqs[-1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
    #dlys = 2*n.pi*n.arange(0, .5/(freqs[1]-freqs[0]) + 1./(freqs[-1]-freqs[0]), 1./(freqs[-1]-freqs[0]))
    sky_maxdly = UMAG / freqs[-1]
    sky_maxbin = n.argmin(n.abs(dlys - sky_maxdly))
    #skymask = n.where(dlys / sky_maxdly < 1., 1., 0)
    skymask = (1 - dlys / (1.*sky_maxdly)).clip(0,1)
    #skymask = (1 - dlys / (2*sky_maxdly)).clip(0,1)
    win = window_hanning(freqs.size)
    #win = window_blackman_harris(freqs.size)
    win /= n.sqrt(n.average(win**2))
    #win = 1.
    zs = 1.420405 / freqs - 1
    #Ls = (freqs * 1e3) * dL_dFQ(zs)
    Ls = (freqs * 1e3) * dL_dFQ(n.average(zs))
    k_pr = dk_du(n.average(zs)) * UMAG
    #k_pl = n.arange(0, .5/(Ls[1]-Ls[0]) + 1./(Ls[-1]-Ls[0]), 1./(Ls[-1]-Ls[0]))
    print 'u:', UMAG, 'FQ:', .5*(FQ1+FQ2), 'z:', n.average(zs)
    print 'L:', Ls[-1]-Ls[0], 'dL:', Ls[1]-Ls[0]
    print '-'*20
    k_pl = 2*n.pi*n.arange(0, .5/(Ls[1]-Ls[0]) + 1./(Ls[-1]-Ls[0]), 1./(Ls[-1]-Ls[0]))
    ks = n.sqrt(k_pl**2 + k_pr**2)
    sky_maxk = ks[sky_maxbin]
    #ks = ks.clip(ks[1]/4, n.Inf)
    #NOISE_LEV = .00015
    NOISE_LEV = .00005
    #NOISE_LEV = .00001
    #NOISE_LEV = .0003
    TPSRC_LEV = .6
    TSYNC_LEV = 1./200

    Tsync = tsky_sync(240., -2.5)
    Tpsrc = tsky_psrc(TPSRC_LEV, skymask, -1)
    Teor_late = pspec2spec(late_eor_pspec(ks))
    Teor_early = pspec2spec(early_eor_pspec(ks))
    Tsky = Tsync * TSYNC_LEV + Tpsrc
    pspec0 = n.abs(n.fft.rfft(Tsky*win))**2 / Tsky.size

    #p.subplot(121)
    #p.semilogy(freqs, n.abs(Tsky), 'k--')
    #p.semilogy(freqs, n.abs(Teor_late), 'k:')
    #p.semilogy(freqs, n.abs(Teor_early), 'k-.')
    #p.xlabel('Frequency (MHz)')
    #p.ylabel('K')
    #p.grid()

    #p.subplot(122)
    _sum1,_wgt1 = 0,0
    _sum2,_wgt2 = 0,0
    _sum3,_wgt3 = 0,0

    for i in range(40):
        Tpsrc = tsky_psrc(TPSRC_LEV, skymask, -1)
        Teor_late = pspec2spec(late_eor_pspec(ks))
        Teor_early = pspec2spec(early_eor_pspec(ks))
        Tsky_none = Tsync * TSYNC_LEV + Tpsrc ; Tsky_none -= .9*n.average(Tsky_none)
        Tsky_late = Tsync * TSYNC_LEV+ Tpsrc + Teor_late ; Tsky_late -= .9*n.average(Tsky_late)
        Tsky_early = Tsync * TSYNC_LEV+ Tpsrc + Teor_early ; Tsky_early -= .9*n.average(Tsky_early)
        noise = NOISE_LEV * n.random.normal(scale=Tsync)
        spec1 = n.fft.rfft((Tsky_none + noise)*win)
        spec2 = n.fft.rfft((Tsky_late + noise)*win)
        spec3 = n.fft.rfft((Tsky_early + noise)*win)
        #spec2 = n.fft.rfft(pspec2spec(late_eor_pspec(ks)))
        _sum1 += n.abs(spec1)**2
        _wgt1 += spec1.size
        _sum2 += n.abs(spec2)**2
        _wgt2 += spec2.size
        _sum3 += n.abs(spec3)**2
        _wgt3 += spec3.size

    pspec1 = _sum1 / _wgt1
    pspec2 = _sum2 / _wgt2
    pspec3 = _sum3 / _wgt3

    _ks,_pspec0 = rebin_log(ks[3:], pspec0[3:])
    _ks,_pspec1 = rebin_log(ks[3:], pspec1[3:])
    _ks,_pspec2 = rebin_log(ks[3:], pspec2[3:])
    _ks,_pspec3 = rebin_log(ks[3:], pspec3[3:])
    _ks = n.concatenate([ks[:3], _ks])
    _pspec0 = n.concatenate([pspec0[:3], _pspec0])
    _pspec1 = n.concatenate([pspec1[:3], _pspec1])
    _pspec2 = n.concatenate([pspec2[:3], _pspec2])
    _pspec3 = n.concatenate([pspec3[:3], _pspec3])

    color = 'kbgrcmy'[cnt%7]
    (fa1,fb1,fc1) = fit_pwrlaw_noise(_ks,_pspec1)
    (fa2,fb2,fc2) = fit_pwrlaw_noise(_ks,_pspec2)

    if True:
        #p.loglog(_ks, _pspec0, 'k--')
        p.loglog(_ks, 1e6*_pspec1, color+'.')
        p.loglog(_ks, 1e6*n.abs(_pspec1-fc1), color+':')
        p.loglog(_ks, 1e6*_pspec2, color+'^')
        p.loglog(_ks, 1e6*n.abs(_pspec2-fc2), color+'-')
        #p.loglog(_ks, _pspec3, 'g')
        p.loglog(ks,1e6*early_eor_pspec(ks), 'k--')
        p.loglog(ks,1e6*late_eor_pspec(ks), 'k--')
        p.loglog([sky_maxk, sky_maxk], [1e0, 1e6], 'k--')
        p.ylabel('$mK^2$')
        p.ylim(1e0,1e6)
    else:
        #p.loglog(_ks, n.sqrt(_pspec0), 'k--')
        p.loglog(_ks, n.sqrt(_pspec1), color+':')
        p.loglog(_ks, n.sqrt(_pspec2), color+'-')
        #p.loglog(_ks, n.sqrt(_pspec3), 'g')
        #p.loglog(ks,n.sqrt(early_eor_pspec(ks)), 'k-.')
        p.loglog(ks,n.sqrt( late_eor_pspec(ks)), 'k:')
        p.ylabel('$K$')
        p.ylim(1e-3,1e0)
  p.title(str(UMAG))
  #p.xlim(1e-2,1e1)
  p.xlim(3e-2,3e0)
  p.xlabel(r'$k (h\ {\rm Mpc})^{-1}$')
  p.grid()
p.show()
    
