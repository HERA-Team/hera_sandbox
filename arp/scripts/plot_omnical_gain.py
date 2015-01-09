#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

#sep = 'sep34'
sep = 'sep6'

times = []
gains = []
for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    times.append(npz['lsts'])
    gains.append(npz['gain_'+sep])
times = n.concatenate(times)
gains = n.concatenate(gains)
p.plot(times, gains, ',')
#p.plot(gains)
p.show()
