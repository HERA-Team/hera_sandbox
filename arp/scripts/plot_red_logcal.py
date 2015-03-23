#! /usr/bin/env python
import numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
o.add_option('-a','--ant', type='int', help='Antenna to plot')
o.add_option('-p','--pol', help='Polarization to plot')
opts,args = o.parse_args(sys.argv[1:])

key1 = '%d,%s' % (opts.ant, opts.pol)
key2 = '%d,%s' % (4, opts.pol)
keys = [key1,key2,'chi2_lin','sep6','iters']
d = {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    for key in keys:
        if not d.has_key(key): d[key] = []
        d[key].append(f[key])

for k in keys: d[k] = n.concatenate(d[k], axis=0)
for k in keys: d[k] = n.where(n.isnan(d[k]), 0, d[k])
d12 = d[key2] * n.conj(d[key1])
d12[:,:30] = 0
d12[:,170:] = 0
w12 = n.where(d12 == 0, 0, 1)
_d12 = C.arp.clean_transform(d12, w=w12, clean=1e-5, window='blackman-harris')
_d12 = n.fft.fftshift(_d12)
p.subplot(141)
#C.arp.waterfall(d12, mode='lin', mx=1.25, drng=.5)
C.arp.waterfall(_d12, mode='log', drng=3)
p.colorbar(shrink=.5)
p.subplot(142)
#C.arp.waterfall(d12, mode='phs', mx=n.pi, drng=2*n.pi)
C.arp.waterfall(d12, mode='phs', mx=.5, drng=1)
p.colorbar(shrink=.5)

p.subplot(143)
C.arp.waterfall(d['chi2_lin'], mode='log')
#C.arp.waterfall(d['iters'], mode='log')
p.colorbar(shrink=.5)

p.subplot(144)
#C.arp.waterfall(d['sep6'], mode='real')
C.arp.waterfall(d['iters'], mode='lin')
p.colorbar(shrink=.5)

p.show()
