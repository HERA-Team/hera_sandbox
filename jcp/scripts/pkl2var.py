#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse, qPickle

o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True)
o.add_option('-P', '--pkl', dest='pkl', help='The pickle to load.')
opts, args = o.parse_args(sys.argv[1:])

pickle = qPickle.load(opts.pkl)

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)

fluxes = []
for src in srclist:
    for flux in pickle['flux'][src]:
        print flux, (10**pickle['flux'][src][flux]).real
        fluxes.append((10**pickle['flux'][src][flux]).real)


print (n.array(fluxes)).var()
