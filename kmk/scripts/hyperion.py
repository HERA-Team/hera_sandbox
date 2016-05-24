#! /usr/bin/env python

import aipy as a
import numpy as np
import pylab as pl

Tsky = 4000 # K
T21 = 20e-3 # K
Trx = 100 # K
Tabs = 300 # K
B = 100 # MHz
t = 300e6 # us, 5 min
# NOTE: changing the scaling doesn't seem to change the Trms??
beam_scaling = 0.10 

im = a.img.Img(200,0.5)
tx,ty,tz = im.get_top()
omega = np.where(tz.mask,0,1)
pit = np.where(tx**2 + ty**2 < beam_scaling,1,0).filled(0)
omegaPrime = omega * pit # baffle-attenuated sky response coverage

T = Trx * omega + Tabs * (omega - omegaPrime) + (Tsky + T21) * omegaPrime
Trms = T / np.sqrt(2*B*t)
print "rms T = %f" % Trms.max()
print "global signal T = %f" % T21

uv_sensitivity = np.fft.fft2(omegaPrime)
# according to Presley 2015, want antennas separated by 1 wavelength
# 1 wavelength == 1 unit in uv-coordinates
# 2 bins in array == 1 unit in uv-coordinates
print uv_sensitivity[2,0] / uv_sensitivity[0,0]

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
