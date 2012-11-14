#!/usr/bin/env python
#
#  radio_lunar_eclipse_find.py
#  
#
#  Created by Danny Jacobs on 12/21/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m,ephem
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('-t',type='float',
    help='Start time (jd)')
o.add_option('--dt',type='float',
    help='Length of search (days)')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])


aa = a.cal.get_aa(opts.cal,n.array([0.1]))
srcs,fluxes,catalogs = a.scripting.parse_srcs(opts.src,opts.cat)
cat = a.cal.get_catalog(opts.cal,srcs,fluxes,catalogs)
M = ephem.Moon()
curt = 0
septest=0.5*a.img.deg2rad#seperation to look for (degrees)
for t in n.linspace(opts.t,opts.t+opts.dt,num=opts.dt*1000):
    for src in cat:
        aa.set_jultime(t)
        cat[src].compute(aa)
        M.compute(aa)
        try: 
            if ephem.separation(cat[src],M)<septest:
                print 'Date = ',aa.date,'\t',
                print 'el = ', cat[src].alt,'\t',
                print 'sep = ', ephem.separation(cat[src],M),'\t',
                print 'name = ', cat[src].src_name,'\t',
                print 'Transit =',cat[src].transit_time,
                print aa.get_jultime(),'\t',
                print "Flux = ",cat[src]._jys
        except(TypeError):continue
        
