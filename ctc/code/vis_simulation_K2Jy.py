#!/usr/bin/env python

"""

NAME:
     vis_simulation_K2Jy.py
PURPOSE:
     Converts UV files with temperature [K] visibilities to [Jy] visibilities
EXAMPLE CALL:
     ./vis_simulation_K2Jy.py --uvfile <path>
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
o.set_usage('vis_simulation_K2Jy.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--uvold', dest='uvold', default='/Users/carinacheng/capo/ctc/tables/vis_simulation_1002.uv', 
             help='Path where UV file is.')
o.add_option('--uvnew', dest='uvnew', default='/Users/carinacheng/capo/ctc/tables/test.uv', 
             help='Path and name of outputted UV file.')
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is .1/203')
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
opts, args = o.parse_args(sys.argv[1:])



def mfunc(uv,p,d):

    d *= 10**3 #[K] to [mK]

    for i in range(len(freqs)):
        d[i] = d[i]/capo.pspec.jy2T(freqs[i])

    return p,d



freqs = numpy.linspace(opts.sfreq,opts.sfreq+opts.sdf*opts.nchan,num=opts.nchan, endpoint=False) #array of frequencies

uvi = aipy.miriad.UV(opts.uvold)
#p,d,f = uvi.read(raw=True)
#print d
#print mfunc(uvi,p,d)[1]

uvo = aipy.miriad.UV(opts.uvnew, status='new')
uvo.init_from_uv(uvi)
uvo.pipe(uvi,mfunc=mfunc)
print opts.uvold, '->', opts.uvnew


