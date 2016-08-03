#!/usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse,sys
import time

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True,src=True)
o.add_option('--jdrng', action='store', 
              help='range of juliandates to use')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))

if opts.src != None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    
jdrange = map(float, opts.jdrng.split('_'))
jdrange = n.arange(jdrange[0], jdrange[1], 1/24./60.)

ras = {}
dec = {}

srcs = n.array(cat.keys())
print srcs



for t in jdrange:
    aa.set_jultime(t)
    cat.compute(aa)
    az,alt= cat.get_crds('top',ncrd=2)
    up = n.where(alt>0)[0]
    above = alt[up]
    around = az[up]
    print srcs[up], " are visible. " 
    p.plot(around,above, 'o')

p.show()

    

