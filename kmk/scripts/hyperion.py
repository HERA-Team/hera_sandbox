#! /usr/bin/env python

import aipy as a
import numpy as np
import pylab as pl

pl.rc('text', usetex=True)
pl.rc('font', family='serif')

N = 1024 # number of pixels on one axis
B = 5 # channel BW, in MHz
t = 0.5e6 # integration time, in us, 5 hr 
beam_scaling = 0.25

im = a.img.Img(N,0.5)
tx,ty,tz = im.get_top()
omega = np.where(tz.mask,0,1./N**2)
pit = np.where(tx**2 + ty**2 < beam_scaling,1,0).filled(0)
omegaPrime = omega * pit # baffle-attenuated sky response coverage
omegaSum = np.sum(omega)
omegaPrimeSum = np.sum(omegaPrime)

def sysTemp(feedBeam, baffleBeam, Tsky=2000, T21=20e-3, Tabs=300, Trx=100):
    '''
        Tsky = 2000 # sky temp, K
        T21 = 20e-3 # 21cm global signal temp fluctuation, K
        Trx = 100 # receiver temp, K
        Tabs = 300 # absorber temp, K
    '''
    T = (Trx * feedBeam) + Tabs * (feedBeam - baffleBeam) + ((Tsky + T21) * baffleBeam)
    return T / baffleBeam

feedBeam = 1
baffleBeam = np.linspace(0.01, feedBeam, 200)
Trx = np.arange(50, 350, 50)
for i in np.arange(len(Trx)):
    T = sysTemp(feedBeam, baffleBeam, Trx=Trx[i])
    #Trms = T / np.sqrt(2*B*t)
    pl.loglog(baffleBeam, T, label=r"T$_{rx}$ = %d K" % Trx[i])

pl.xlabel(r"Fractional Sky Coverage ($\Omega'/\Omega$)", fontsize=18)
pl.ylabel(r'System Temperature (K)', fontsize=18)
pl.tick_params(axis='both', direction='in', color='black', labelsize=18)
pl.grid(b=True, which='both')
pl.legend()
pl.savefig('sysTemp.png')
pl.show()

#print "rms T = %f" % Trms
#print "global signal T = %f" % T21

#uv_sensitivity = np.fft.fft2(omegaPrime)
#uv_sensitivity = uv_sensitivity / np.amax(uv_sensitivity) # normalize
# according to Presley 2015, want antennas separated by 1 wavelength
# 1 wavelength == 1 unit in uv-coordinates
# 2 bins in array == 1 unit in uv-coordinates
#print uv_sensitivity[2,0] 

#pl.figure(0)
#pl.subplot(121)
#pl.imshow(np.fft.fftshift(omega), extent=(-1,1,-1,1), interpolation="nearest")
#pl.colorbar()
#assert (np.all(omegaPrime == -1*omega))# "these are not the same"
#pl.imshow(np.fft.fftshift(omegaPrime), extent=(-1,1,-1,1), interpolation="nearest", alpha=0.5)
#pl.subplot(122)
#pl.imshow(np.fft.fftshift(np.abs(uv_sensitivity)), extent=(-100,100,-100,100), interpolation="nearest")
#pl.colorbar()
#pl.show()
