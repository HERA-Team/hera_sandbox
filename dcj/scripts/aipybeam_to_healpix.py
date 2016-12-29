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
o.add_option('-f','--freq',dest='freq',default=0.15,type=float,
    help='desired frequency in GHz')
o.add_option('--pol',default='x',
    help='Polarization default=x')
o.add_option('--nside',default=512,
    help='nside at which to perform integral')
o.add_option('--filename',default='beam.fits',
    help="beam file to output default=beam.fits")
opts, args = o.parse_args(sys.argv[1:])

npix = hp.nside2npix(opts.nside)
pix = np.arange(npix)
theta,phi = hp.pix2ang(opts.nside,pix)
#only keep points above the horizon
# phi[theta<90] = phi[theta<90]
# theta[theta<90] = theta[theta<90]
intpix = hp.ang2pix(opts.nside,theta,phi)
x,y,z = hp.pix2vec(opts.nside,intpix)
pixarea = hp.nside2pixarea(opts.nside) #pixel area in radians

print "loading beam model"
aa = a.cal.get_aa(opts.cal, np.array([opts.freq]))

print("evaluating the beam model for %i points"%npix)
A = aa[0].bm_response((x,y,z),pol=opts.pol)
hp.write_map(opts.filename,A)
