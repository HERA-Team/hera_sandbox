import aipy as a, numpy as np
import sys, optparse
from pylab import *
import healpy as hp

"""
Integrate the beam model stored in an aipy cal file
returns Omega, 1/Omega
and the sqrt(1/Omega) which is the "dish_size_in_lambda" expected by 21cmsense
"""

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-f','--freq',dest='freq',default='0.15_0.16',type=str,
    help='frequency range in GHz over which to integrate beam model')
o.add_option('--pol',default='x',
    help='Polarization default=x')
o.add_option('--nside',default=512,
    help='nside at which to perform integral')
opts, args = o.parse_args(sys.argv[1:])

fstart,fstop = map(float,opts.freq.split('_'))
freqs = np.linspace(fstart,fstop)
B = freqs.max()-freqs.min()
df = np.diff(freqs)[0]
freq = np.mean(freqs)
print "loading beam model"
aa = a.cal.get_aa(opts.cal, freqs)
npix = hp.nside2npix(opts.nside)
pix = np.arange(npix)
theta,phi = hp.pix2ang(opts.nside,pix)
#only keep points above the horizon
phi[theta<90] = phi[theta<90]
theta[theta<90] = theta[theta<90]
intpix = hp.ang2pix(opts.nside,theta,phi)
x,y,z = hp.pix2vec(opts.nside,intpix)
pixarea = hp.nside2pixarea(opts.nside) #pixel area in radians


#nomenclature after Parsons 2016 (1503.05564v2)
print "integrating over bandwidth"
A = aa[0].bm_response((x,y,z),pol=opts.pol)**2  #A is defined as the power beam, there is one term of it in the visibility equation
Omega_pp = np.sum(A**2) * pixarea * df / B
print "Omega_pp = ",np.round(Omega_pp,2)
Omega_p = np.sum(A)*pixarea * df / B
print "Omega_p = ",Omega_p
print "Omega_p^2/Omega_pp = ",Omega_p**2/Omega_pp 


print("computing the beam at %s GHz"%freq)
aa[0].beam.afreqs = np.array([freq])
A = aa[0].bm_response((x,y,z),pol=opts.pol)**2
Omega = np.sum(A) * pixarea
print "beam_size_in_wavelengths",1/np.sqrt(Omega)
wl = 3e8/(freq*1e9)
print "Effective Area",np.round(wl**2/Omega,2),"m"
