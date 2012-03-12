#!/usr/bin/env python
#
#  send_sun_to_topcat.py
#  
#
#  Created by Danny Jacobs on 4/28/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

t = float(args[0])
aa = a.cal.get_aa(opts.cal,0.0001,0.1,1)
sun = ephem.Sun()
aa.set_jultime(t)
sun.compute(aa)


VO = atpy.Table()
VO.add_column('Ra',[repr(sun.ra)],dtype='float',units='radians')
VO.add_column('Dec',[repr(sun.dec)],dtype='float',units='radians')
