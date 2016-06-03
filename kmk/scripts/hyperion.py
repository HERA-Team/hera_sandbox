#! /usr/bin/env python

import aipy as a
import numpy as np
import pylab as pl

N = 200 # number of pixels on one axis
Tsky = 4000 # K
T21 = 20e-3 # K
Trx = 100 # K
Tabs = 300 # K
B = 5 # channel BW, in MHz
t = 5*3600e6 # integration time, in us, 5 hr 
beam_scaling = 0.25

im = a.img.Img(N,0.5)
tx,ty,tz = im.get_top()
omega = np.where(tz.mask,0,1)
pit = np.where(tx**2 + ty**2 < beam_scaling,1,0).filled(0)
omegaPrime = omega * pit # baffle-attenuated sky response coverage

T = np.sum(Trx * omega) + np.sum(Tabs * (omega - omegaPrime)) + np.sum((Tsky + T21) * omegaPrime)
T = T / N**2
Trms = T / np.sqrt(2*B*t)
print "rms T = %f" % Trms
print "global signal T = %f" % T21

uv_sensitivity = np.fft.fft2(omegaPrime)
uv_sensitivity = uv_sensitivity / np.amax(uv_sensitivity) # normalize
# according to Presley 2015, want antennas separated by 1 wavelength
# 1 wavelength == 1 unit in uv-coordinates
# 2 bins in array == 1 unit in uv-coordinates
print uv_sensitivity[2,0] 

pl.figure(0)
pl.subplot(121)
pl.imshow(np.fft.fftshift(omega), extent=(-1,1,-1,1), interpolation="nearest")
pl.colorbar()
#assert (np.all(omegaPrime == -1*omega))# "these are not the same"
pl.imshow(np.fft.fftshift(omegaPrime), extent=(-1,1,-1,1), interpolation="nearest", alpha=0.5)
pl.subplot(122)
pl.imshow(np.fft.fftshift(np.abs(uv_sensitivity)), extent=(-100,100,-100,100), interpolation="nearest")
pl.colorbar()
pl.show()
