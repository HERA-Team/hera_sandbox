#! /usr/bin/env python
import aipy as a
import numpy as np
import matplotlib.pylab as pl

# NOTE: optimal baseline length formula given in Presley et al 2015 eq. 9

freq = np.array([.100]) # in GHz
c = 3e10 # in cm/s
wvlen = c / freq # observed wavelength
b = 4e2 / a.const.len_ns # baseline between antennae, in nanoseconds
b = np.array([b, 0, 0])

# fill h with flat temperature across whole sky (DC signal)
h = a.healpix.HealpixMap(nside=64)
h.map = np.ones(shape = h.map.shape)
h.map = h.map*4*np.pi/h.npix()
x, y, z = h.px2crd(np.arange(h.npix())) #topocentric

# primary beam of antenna
beam = a.amp.Beam2DGaussian(freqs = freq, xwidth=np.pi, ywidth=np.pi)
beam = beam.response(xyz=(x,y,z))[0,:]

# attenuate sky signal by primary beam
obs = beam * h.map

phs = np.exp(-2j*np.pi*freq*(b[0]*x + b[1]*y + b[2]*z))

vis = np.sum(np.where(z>0,obs*phs,0))

print vis


#ants = []
#ants.append(a.phs.Antenna(0,0,0, beam, delay=0, offset=0))
#ants.append(a.phs.Antenna(0,b,0, beam, delay=0, offset=0))

#aa = a.ant.AntennaArray(ants=ants, location=("18:20:39", "-66:45:10"))
