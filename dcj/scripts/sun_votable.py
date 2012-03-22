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
import atpy

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

tmin = float(args[0])
tmax = float(args[1])
aa = a.cal.get_aa(opts.cal,0.0001,0.1,1)
sun = ephem.Sun()

pos = []
for t in n.linspace(tmin,tmax):
    aa.set_jultime(t)
    sun.compute(aa)
    pos.append((t,float(sun.ra),float(sun.dec)))

VO = atpy.Table()
VO.add_column('Ra',[l[1] for l in pos],dtype='float',unit='radians')
VO.add_column('Dec',[l[2] for l in pos],dtype='float',unit='radians')
VO.add_column('_Ra',[l[1]*180/n.pi for l in pos],dtype='float',unit='deg')
VO.add_column('_Dec',[l[2]*180/n.pi for l in pos],dtype='float',unit='deg')

VO.add_column('time',[l[0] for l in pos],dtype='float',unit='jd')
VO.write('sun_'+'-'.join(args)+'.vot')
