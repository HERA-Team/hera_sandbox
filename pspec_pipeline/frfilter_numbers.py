#!/usr/bin/env python
import numpy as n
import capo as C
import aipy as a
from matplotlib import pylab as p
import scipy
import capo.frf_conv as fringe


def beam_area(hmap):
    return 4*n.pi*n.sum(hmap)/h.npix()

#get antenna array
aa = a.cal.get_aa('psa6622_v003', n.array([.159])) #XXX hard-coded

h = a.healpix.HealpixMap(nside=64) #healpix map for the beam
xyz = h.px2crd(n.arange( h.npix() ), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz) #rotate the coordinated system to be centered on the array. This is equatorial centered at the array.
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmI = 0.5 * (_bmx**2 + _bmy**2)
bmI = n.where(tz > 0, bmI, 0) # only use beam values above the horizon.

bl = aa.get_baseline(0,26,'r') * .151 #XXX hard-coded baseline length in frequency.
print aa.get_baseline(0,26,'r')
fng = fringe.mk_fng(bl,xyz)

#get the fringe rate filter in frf_conv. aa only has one channel in it.
frp, bins = fringe.aa_to_fr_profile(aa, (1,4), 0) #XXX hard-coded
tbins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[0],alietal=False)
#frp = fringe.fir_to_frp(firs)

## Old way of calculating FRP,FIRS,etc. ##
#frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,26)) 
#bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,26), ref_chan=0) 
#tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)

#need to renormalize to get the proper scaling. firs are properly normalized.
#firs = firs[0]
#frfs = n.fft.fftshift(n.fft.fft(n.fft.ifftshift(firs), axis=-1))
#get weights.
wgts = scipy.interpolate.interp1d(bins, frp, kind='linear')
fng_wgt = wgts(fng) #gets weightings at fringes on the sky.
fng_bm = bmI * fng_wgt
#flat weighting determined by the maximum possible fringe rate for a baseline 
#and 0.
skypass = n.where(n.abs(bins) < n.max(n.abs(fng)), 1., 0)
# This is the noise level after being filtered and averaged by the application of this filter.
crit_noise_lev = n.sum(skypass**2) / bins.size
#ditto for the optimal weights.
fng_noise_lev = n.sum(n.abs(frp)**2)/ bins.size

print 'integration times PSA128', 31.6*(1/crit_noise_lev), 31.6*(1/fng_noise_lev)
print 'integration times PSA64', 43*(1/crit_noise_lev), 43*(1/fng_noise_lev)
print 'beams (flat) (optimal) [power]', beam_area(bmI), beam_area(fng_bm)
print 'beams (flat) (optimal) [power^2]', beam_area(bmI**2), beam_area(fng_bm**2)
#signal loss is the ratio of the beams.
print 'signal loss', n.sum(fng_bm**2)/n.sum(bmI**2)
#noise reduction is the ratio of the standard deviations
print 'noise reduction', crit_noise_lev/fng_noise_lev
print 'sensitivity', (crit_noise_lev/fng_noise_lev)**.5 / (n.sum(bmI**2)/n.sum(fng_bm**2))

