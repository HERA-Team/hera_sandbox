#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, optparse

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

avg = n.mean(a.img.from_fits(args[0])[0])
for map in args:
    map = a.img.from_fits(map)[0]
    gmean *= map*1e8

#gmean = gmean**(1./len(gmeans))

print gmean

a.img.to_fits('gmean.fits',gmean)
