#!/usr/bin/env python

"""

NAME: 
      pspec2vis_sim_wrapper.py
PURPOSE:
      -Wrapper code for making GSM and going from a power spectrum to visibility
EXAMPLE CALL:
      ./pspec2vis_sim_wrapper.py --nchan 203 
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import optparse
import os, sys


#options

o = optparse.OptionParser()
o.set_usage('vis_simulation_v2.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--map', dest='map', default='make',
             help='Directory where GSM files are (labeled gsm1001.fits, gsm1002.fits, etc.). Include final / when typing path. Default makes GSM files.')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/test.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
o.add_option('--inttime', dest='inttime', default=10., type='float',
             help='Integration time (s). Default is 10.') 
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is .1/203')
o.add_option('--startjd', dest='startjd', default=2454500., type='float',
             help='Julian Date to start observation.  Default is 2454500')
o.add_option('--endjd', dest='endjd', default=2454501., type='float',
             help='Julian Date to end observation.  Default is 2454501')
o.add_option('--psa', dest='psa', default='psa898_v003', 
             help='Name of PSA file.')
o.add_option('--lmax', dest='lmax', default=5, type='int',
             help='Maximum l value. Default is 5.')
opts, args = o.parse_args(sys.argv[1:])

#make GSM if needed

if opts.map == 'make':

    print 'MAKING GSM...'

    
    
