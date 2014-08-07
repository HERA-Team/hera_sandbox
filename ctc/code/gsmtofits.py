#!/usr/bin/env python

"""

NAME: 
      gsmtofits.py
PURPOSE: 
      Converts a folder of GSM files to fits files. Need a file called 'filenames' that lists all the .dat files. This file needs to be in the same directory as the .dat files. See description of vis_simulation.py for instructions on making GSM .dat files.
EXAMPLE CALL:
      ./gsmtofits.py --map /Users/carinacheng/capo/ctc/images/gsm/gsm256/
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import asciidata
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('gsmtofits.py [options]')
o.set_description(__doc__)
o.add_option('--map', dest='map', default='/Users/carinacheng/capo/ctc/images/gsm/gsm256/',
             help='Directory where filenames and GSM files are (labeled gsm1001.dat, gsm1002.dat, etc.). Include final / when typing path.')
ops, args = o.parse_args(sys.argv[1:])


root = ops.map
f = asciidata.open(ops.map + 'filenames')
names = f[0].tonumpy()

for ii in range(len(names)):
    
    d = numpy.loadtxt(root+names[ii])

    h = aipy.healpix.HealpixMap(nside=512)
    h.map = d

    h.to_fits((root+names[ii]).replace('.dat','.fits'))

    print names[ii] + ' completed'
             

                
