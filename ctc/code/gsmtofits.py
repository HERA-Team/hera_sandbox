## NAME: 
#         gsmtofits.py
## PURPOSE: 
#         Converts a folder of gsm files to fits files

#------------------------------------------------------------------------------------------------

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import asciidata

root = '/Users/carinacheng/Desktop/Carina/UCBResearch/images/gsm/gsm256/'
f = asciidata.open(root + 'filenames')
names = f[0].tonumpy()

for ii in range(len(names)):
    
    d = numpy.loadtxt(root+names[ii])

    h = aipy.healpix.HealpixMap(nside=512)
    h.map = d

    h.to_fits((root+names[ii]).replace('.dat','.fits'))

    print names[ii] + ' completed'
             

                
