#! /usr/bin/env python
import numpy as n, pylab as p
import sys

for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    for k in npz.files:
        if k.startswith('t'): continue
        p.plot(npz['t'+k], npz[k], ',')
p.show()
