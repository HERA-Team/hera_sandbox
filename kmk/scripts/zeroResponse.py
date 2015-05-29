#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl

# NOTE: optimal baseline length formula given in Presley et al 2015 eq. 9
# fill h with flat temperature across whole sky (DC signal)
h = a.healpix.HealpixMap(nside=64)
h.map = np.ones(shape = h.map.shape)
h.map = h.map*4*np.pi/h.npix()
#temps = np.ones(3)
temps = np.array([0.5, 0.6, 0.8]) # sky temperature
x, y, z = xyz = h.px2crd(np.arange(h.npix())) #topocentric

freq = np.array([.120, .150, .180]) # in GHz
c = 3e10 # in cm/s
wvlen = c / freq # observed wavelength

# import array parameters
aa = a.cal.get_aa('psa6622_v001', freq)

beam = aa[0].bm_response(xyz, pol='x')**2
beam = beam[0]

bx, by, bz = aa.get_baseline(0, 103)
print bx, by, bz

vis = np.zeros(shape = temps.shape, dtype=np.complex) 

# attenuate sky signal by primary beam
for i in range(len(temps)):
    obs = temps[i] * beam * h.map
    phs = np.exp(-2j*np.pi*freq[i]*(bx*x + by*y + bz*z))
    vis[i] = np.sum(np.where(z>0,obs*phs,0))
    i += 1

print vis
