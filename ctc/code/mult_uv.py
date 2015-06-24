#!/usr/bin/env python

"""

NAME:
     mult_uv.py
PURPOSE:
     Multiplies UV file data by some factor
EXAMPLE CALL:
     ./mult_uv.py --uvold <path> --uvnew <path> --factor 100
AUTHOR:
     Carina Cheng

"""

import pylab
import numpy
import optparse
import pyfits
import aipy
import capo
import matplotlib.pyplot as plt
import os, sys

o = optparse.OptionParser()
o.set_usage('mult_uv.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--uvold', dest='uvold', default='/Users/carinacheng/capo/ctc/tables/vis_simulation_1002.uv', 
             help='Path where UV file is.')
o.add_option('--uvnew', dest='uvnew', default='/Users/carinacheng/capo/ctc/tables/test.uv', 
             help='Path and name of outputted UV file.')
o.add_option('--factor',dest='factor', default=1000, type='float',
             help='Factor to multiply UV file data by. Default is 1000.')
opts, args = o.parse_args(sys.argv[1:])


def mfunc(uv,p,d):

    d *= opts.factor

    return p,d

uvi = aipy.miriad.UV(opts.uvold)
#p,d,f = uvi.read(raw=True)
#print d
#print mfunc(uvi,p,d)[1]

uvo = aipy.miriad.UV(opts.uvnew, status='new')
uvo.init_from_uv(uvi)
uvo.pipe(uvi,mfunc=mfunc)
print opts.uvold, '->', opts.uvnew


