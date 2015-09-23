#! /usr/bin/env python
'''
Implements a Fringe Rate Filter useing and FIR filter.
'''
import aipy as a
import capo
import numpy as n, pylab as p
from numpy.fft import ifft, fftshift, ifftshift, fftfreq, fft
import scipy.interpolate
import optparse,sys

DEFAULT_FRBINS = n.arange(-.01+5e-5/2,.01,5e-5) # Hz
DEFAULT_WGT = lambda bm: bm**2
DEFAULT_IWGT = lambda h: n.sqrt(h)

def mk_fng(bl, eq):
    '''Return fringe rates given eq coords and a baseline vector (measured in wavelengths) in eq coords'''
    return -2*n.pi/a.const.sidereal_day * n.dot(n.cross(n.array([0,0,1.]),bl), eq)

#fringe used in ali et.al to degrade optimal fringe rate filter.
#def mk_fng_alietal(bl, ex, ey, ez):
#    return 2*n.pi/a.const.sidereal_day * (bl[0]*ex + bl[1]*ey * n.sqrt(1-ez**2))

def fr_profile(bm, fng, bins=DEFAULT_FRBINS, wgt=DEFAULT_WGT, iwgt=DEFAULT_IWGT):
    '''Return the fringe-rate profiel (binning the beam by fringe rate).'''
    h, _ = n.histogram(fng, bins=bins, weights=wgt(bm)) 
    h = iwgt(h)
    h /= h.max()
    #bins given to histogram are bin edges. Want bin centers.
    bins = 0.5 * (bins[:-1] + bins[1:]) 
    return h, bins

def gauss(cenwid, bins): return n.exp(-(bins-cenwid[0])**2/(2*cenwid[1]**2))
def tanh(x, p, w, C = 1.0, a=1.0): return (C/2.) * (1 + a*n.tanh( (x-p)/(2*w)))
def mdl_wrap(prms, frp, bins, maxfr, mdl): return n.sum((frp - n.where(bins > maxfr,0,mdl(prms,bins)))**2)

def fit_mdl(frp, bins, maxfr, mdl=gauss, maxfun=1000, ftol=1e-6, xtol=1e-6, startprms=(.001,.0001), verbose=False):
    '''Fit a parametrized model to the fringe-rate profile frp.  h is weights for each fringe rate in bins,
    maxfr is the maximum fringe rate (which is treated independently of bins).'''
    bestargs, score = a.optimize.fmin(mdl_wrap, x0=startprms, args=(frp,bins,maxfr,mdl),
        full_output=1, disp=0, maxfun=maxfun, maxiter=n.Inf, ftol=ftol, xtol=xtol)[:2]
    if verbose: print 'Final prms:', bestargs, 'Score:', score
    return bestargs

# XXX wgt and iwgt seem icky
def hmap_to_fr_profile(bm_hmap, bl, lat, bins=DEFAULT_FRBINS, wgt=DEFAULT_WGT, iwgt=DEFAULT_IWGT):
    '''For a healpix map of the beam (in topocentric coords), a bl (in wavelengths, eq coords), 
    and a latitude (in radians), return the fringe-rate profile.'''
    eq = bm_hmap.px2crd(n.arange(bm_hmap.npix()), ncrd=3) # equatorial coordinates
    eq2zen = a.coord.eq2top_m(0., lat)
    top = n.dot(eq2zen, eq)
    bm = bm_hmap[(top[0], top[1], top[2])]
    fng = mk_fng(bl,eq)
    return fr_profile(bm, fng, bins=bins, wgt=wgt, iwgt=iwgt)
    
def aa_to_fr_profile(aa, (i,j), ch, pol='I', bins=DEFAULT_FRBINS, wgt=DEFAULT_WGT, iwgt=DEFAULT_IWGT, nside=64):
    '''For an AntennaArray, for a baseline indexed by i,j, at frequency fq, return the fringe-rate profile.'''
    fq = aa.get_freqs()[ch]
    h = a.healpix.HealpixMap(nside=nside)
    eq = h.px2crd(n.arange(h.npix()), ncrd=3)
    top = n.dot(aa._eq2zen, eq)
    fng = mk_fng(aa.get_baseline(i,j,'r')*fq, eq)
    # XXX computing bm at all freqs, but only taking one
    _bmx = aa[0].bm_response((top), pol='x')[ch]; _bmx = n.where(top[2] > 0, _bmx, 0)
    _bmy = aa[0].bm_response((top), pol='y')[ch]; _bmy = n.where(top[2] > 0, _bmy, 0)
    if   pol == 'xx': bm = _bmx * _bmx.conj()
    elif pol == 'yy': bm = _bmy * _bmy.conj()
    elif pol == 'xy': bm = _bmx * _bmy.conj()
    elif pol == 'yx': bm = _bmy * _bmx.conj()
    elif pol ==  'I': bm = .5 * (_bmx*_bmx.conj() + _bmy*_bmy.conj())
    elif pol ==  'Q': bm = .5 * (_bmx*_bmx.conj() - _bmy*_bmy.conj())
    elif pol ==  'U': bm = .5 * (_bmx*_bmy.conj() + _bmy*_bmx.conj())
    elif pol ==  'V': bm = .5 * (_bmx*_bmy.conj() - _bmy*_bmx.conj())
    return fr_profile(bm, fng, bins=bins, wgt=wgt, iwgt=iwgt)

#wgt = scipy.interpolate.interp1d(bins, h_I, kind='linear')
#max_fr = n.sqrt(n.dot(bl,bl)) * freqs[-1]*2*n.pi / a.const.sidereal_day
#fr_bins = n.arange(-.5/42.8, .5/42.8-fr_bin_res, fr_bin_res) #in Hz!
#fr_bins = n.fft.fftshift(n.fft.fftfreq(nbins, inttime))

# XXX write a function that generates bins from inttime and time window for fir

def fir_to_frp(fir,tbins=None):
    '''Transform a fir (time domain fr filter) to a fringe rate profile.
       fir: array of fringe rate profile. 
       tbins: Corresponding time bins of filter. If None, doesnt return ffringe rates.
    '''
    fir = ifftshift(fir, axes=-1)
    frp = fft(fir, axis=-1)
    frp = fftshift(frp, axes=-1)
    if tbins is not None: return frp, fftshift(fftfreq(tbins.size, tbins[1]-tbins[0]))
    else: return frp

def frp_to_fir(frp, fbins=None):
    '''Transform a fringe rate profile to a fir filter.'''
    frp = ifftshift(frp,axes=-1)
    fir = ifft(frp, axis=-1)
    fir = fftshift(fir, axes=-1)
    if fbins is not None: return fir, fftshift(fftfreq(fbins.size, fbins[1] - fbins[0]))
    else: return fir
    

def frp_to_firs(frp0, bins, fqs, fq0=.150, limit_maxfr=True, limit_xtalk=True, fr_xtalk=.00035, maxfr=None,
        mdl=gauss, maxfun=1000, ftol=1e-6, xtol=1e-6, startprms=(.001,.0001), window='blackman-harris', verbose=False):
    ''' Take a fringe rate profile at one frequency, fit an analytic function and extend 
        to other frequencies. 
        frp0: fringe rate profile at a single frequency. 
        bins: fr bins that correspind to frp0.
        fqs: Frequencies to extend fiter to. 
        fq0: Frequency at which frp0 is made for.
        limit_maxfr: cut of fringe rates above maximum possible fringe rate.
        fr_xtalk: Threshold for removing crosstalk. 
        mdl: a function to fit the fringe rate profile too. gaussian for default.
    '''
    if maxfr is None: maxfr = bins[n.argwhere(frp0 != 0).max()] # XXX check this
    prms0 = fit_mdl(frp0, bins, maxfr, mdl=mdl,maxfun=maxfun,ftol=ftol,xtol=xtol,startprms=startprms,verbose=verbose)
    prms0 = n.array(prms0)
    if limit_maxfr:
        def limit_maxfr(fq): return tanh(bins,maxfr/fq0*fq,1e-5,a=-1.)
    else:
        def limit_maxfr(fq): return 1
    if limit_xtalk: limit_xtalk = tanh(bins,fr_xtalk,1e-5,a=1.)
    else: limit_xtalk = 1
    frps = n.array([mdl(prms0*fq/fq0,bins) * limit_maxfr(fq) * limit_xtalk for i,fq in enumerate(fqs)])
    tbins = fftshift(fftfreq(bins.size, bins[1]-bins[0]))
    firs = frp_to_fir(frps) 
    #frps = ifftshift(frps, axes=-1)
    #firs = ifft(frps, axis=-1)
    #firs = fftshift(firs, axes=-1)
    firs *= a.dsp.gen_window(bins.size, window)
    firs /= n.sum(n.abs(firs),axis=1).reshape(-1,1) # normalize so that n.sum(abs(fir)) = 1
    return tbins, firs
