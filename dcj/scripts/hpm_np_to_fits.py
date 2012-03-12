#!/usr/bin/env python
#
#  hpm_np_to_fits.py
#  
#
#  Created by Danny Jacobs on 3/19/10.
#  PAPER Project
#
"""
convert a txt file with a list of healpix values into a healpix fits format

"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,healpy as hp

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--nest',dest='nest',action='store_true',
    help="Switch from default ring mode to nest.")
opts, args = o.parse_args(sys.argv[1:])

for file in args:
    outfile = '.'.join(file.split('.')[:-1])+'.fits'
    hpm = n.loadtxt(file,comments='#')
    print file+' > '+outfile
    hp.write_map(outfile,hpm,nest=opts.nest)