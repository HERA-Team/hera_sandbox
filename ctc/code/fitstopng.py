#!/usr/bin/env python

"""

NAME: 
      fitstopng.py
PURPOSE: 
      Converts a directory of pspec .fits images to .png images
EXAMPLE CALL:
      ./fitstopng.py --pspec /Users/carinacheng/capo/ctc/images/pspecs/pspec50lmax200
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('fitstopng.py [options]')
o.set_description(__doc__)
o.add_option('--pspec', dest='pspec', default='/Users/carinacheng/capo/ctc/images/pspecs/pspec50lmax200',
             help='Directory where pspec images are contained.')
opts, args = o.parse_args(sys.argv[1:])

#path = '/Users/carinacheng/capo/ctc/images/pspecs/' + opts.pspec
path = opts.pspec

for root, dirs, filenames in os.walk(path):
    for f in filenames:
        if f[9:] == '.fits':
            newname = f.replace('.fits','.png')
            command = 'plot_map.py -m real -o ' + str(newname) + ' ' + str(f)
            os.chdir(path)
            os.system(command)


