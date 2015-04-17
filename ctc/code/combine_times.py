#! /usr/bin/env python

"""

NAME:
    combine_times.py.
PURPOSE:
    Combines a bunch of UV files of continuous time chunks into one UV file
EXAMPLE CALL:
    ./combine_times.py --uvnew <path>
AUTHOR:
    Carina Cheng

"""

import aipy
import numpy
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('combine_times.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--uvnew', default='/Users/carinacheng/capo/ctc/tables/test.uv',
             help='Path and name of outputted UV file.')
opts,args = o.parse_args(sys.argv[1:])

names=[]

for uvfile in args:

    names.append(uvfile)
    #uvi = aipy.miriad.UV(uvfile)
    #print uvfile,'->',opts.uvnew

print names
