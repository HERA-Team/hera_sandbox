#!/usr/bin/env python
#
#  find_facet.py
#  
#
#  Created by Danny Jacobs on 1/19/10.
#  PAPER Project
#
"""
Finds facets for the input sources.
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem

o = optparse.OptionParser()
o.set_usage('find_facet.py -s srcs facets*.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True)
o.add_option('--sep',dest='sep',default=5,type='float',
   help='Selection criteria for source distance from phase center. [5 degrees]')
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')
opts, args = o.parse_args(sys.argv[1:])
opts.sep *= a.img.deg2rad

def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)

srcs,coff,catalogs = a.scripting.parse_srcs(opts.src,opts.cat)
cat = a.src.get_catalog(srcs=srcs,catalogs=catalogs)
update_pos(cat)
for file in args:
    data,kwds = a.img.from_fits(file)
    center = a.phs.RadioFixedBody(kwds['ra']*n.pi/180,
            kwds['dec']*n.pi/180)
    ephem.FixedBody.compute(center,ephem.J2000)

    for name,src in cat.iteritems():
        cursep = ephem.separation(src,center)
        if cursep<opts.sep:
            if opts.verb: print src.src_name,
            if opts.verb: print "%2.2f"%(cursep/a.img.deg2rad,),
            print file
            continue
    