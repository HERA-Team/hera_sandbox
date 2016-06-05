import healpy as hpy
from astropy.io import fits
import sys,os
import numpy as n,healpy as hpy
"""
load a bunch of healpix fits files from rts and produce a single averaged healpix file
"""
nside=1024

files = sys.argv[1:]
D = n.zeros(hpy.nside2npix(nside))
W = n.zeros(hpy.nside2npix(nside))
print D.shape
for filename in files:
    H = fits.open(filename)
    D[H[1].data['HEALPIX_PIXNUM']] += H[1].data['WEIGHTED PP']
    W[H[1].data['HEALPIX_PIXNUM']] += H[1].data['WEIGHT']
D[W>0] /= W[W>0]
print "len(D)",len(D),12*nside**2
outfile = os.path.commonprefix(files)+'.fits'
print "--->",outfile
hpy.write_map(outfile,D,fits_IDL=False)

