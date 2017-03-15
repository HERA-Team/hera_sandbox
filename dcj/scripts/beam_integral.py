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
o.add_option('-f','--freq',dest='freq',default='0.15',type=str,
    help='frequency range in GHz over which to integrate beam model')
o.add_option('--pol',default='x',
    help='Polarization default=x')
o.add_option('--nside',default=512,type=int,
    help='nside at which to perform integral')
o.add_option('--dB',action='store_true',help='input beam is dB, else assumed to be voltage')
opts, args = o.parse_args(sys.argv[1:])

# fstart,fstop = map(float,opts.freq.split('_'))
# freqs = np.linspace(fstart,fstop)
# B = freqs.max()-freqs.min()
# df = np.diff(freqs)[0]
freq = opts.freq
freqs = np.array([0.15])
print "loading beam model"
if len(args)==1:
    print "loading healpix fits file",args[0]
    beammap = hp.read_map(args[0])
    if opts.dB:
        beammap = 10**(beammap/20)
    print "using nside of fits file"
    opts.nside = hp.npix2nside(len(beammap))
npix = hp.nside2npix(opts.nside)
pix = np.arange(npix)
theta,phi = hp.pix2ang(opts.nside,pix)
#only keep points above the horizon
phi = phi[theta<np.pi/2]
theta = theta[theta<np.pi/2]
intpix = hp.ang2pix(opts.nside,theta,phi)
x,y,z = hp.pix2vec(opts.nside,intpix)
print z.max(),z.min()

pixarea = hp.nside2pixarea(opts.nside) #pixel area in radians
print "pixel area in steradians",pixarea
#V is defined as the voltage beam, there are two terms of it in the visibility equation
if len(args)==1:
    theta,phi = hp.vec2ang(np.vstack((x,y,z)))
    V = beammap[intpix]/np.max(beammap)

else:
    aa = a.cal.get_aa(opts.cal, freqs)
    V = aa[0].bm_response((x,y,z),pol=opts.pol)
    V /= V.max()
#nomenclature after Parsons 2016 (1503.05564v2)
print V.max(),V.min()
Omega_V = np.sum(V) * pixarea

print "Omega voltage = ",np.round(Omega_V,2)
A = V**2  #A is defined as the power beam, there is one term of it in the visibility equation
Omega_p = np.sum(A)*pixarea
print "Omega_p = ",Omega_p

Omega_pp = np.sum(A**2) * pixarea
print "Omega_pp = ",np.round(Omega_pp,2)
print "Omega_p^2/Omega_pp = ",Omega_p**2/Omega_pp

print "beam size in wavelengths = ", 1/np.sqrt(Omega_p)
