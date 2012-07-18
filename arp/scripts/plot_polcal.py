#! /usr/bin/env python
import numpy as n, pylab as p
import sys, cPickle

taus = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    f = open(filename)
    _taus = cPickle.load(f)
    for i,tau in _taus.items(): taus[i] = taus.get(i,[]) + [tau]
    f.close()

for i,tau in taus.items():
    d = n.array(tau)
    d_avg = n.median(d)
    print i, 'Dly (avg):', d_avg
    d -= d_avg
    
    p.plot(d, label='%d'%i)

p.legend()
p.show()
