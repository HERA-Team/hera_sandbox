#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

npzfile = args[0]
dat = n.load(npzfile)
keys = dat.files[:]
keys.remove('kpl')

for umag in keys:
    p.loglog(n.abs(dat['kpl']),n.abs(n.real(dat[umag])),label=str(umag))

p.legend()
p.show()

