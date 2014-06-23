#! /usr/bin/env python
import numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
o.add_option('-a','--ant', type='int', help='Antenna to plot')
o.add_option('-p','--pol', help='Polarization to plot')
opts,args = o.parse_args(sys.argv[1:])

key1 = '%d,%s' % (opts.ant, opts.pol)
key2 = '%d,%s' % (41, opts.pol)
keys = [key1,key2]
d = {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    for key in keys:
        if not d.has_key(key): d[key] = []
        d[key].append(f[key])

for k in keys:
    d[k] = n.concatenate(d[k], axis=0)
d12 = d[key1] #* n.conj(d[key1])
p.subplot(121)
C.arp.waterfall(d12, mode='log')
p.colorbar(shrink=.5)
p.subplot(122)
C.arp.waterfall(d12, mode='phs', mx=n.pi, drng=2*n.pi)
p.colorbar(shrink=.5)

p.show()
