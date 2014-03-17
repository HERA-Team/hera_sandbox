#! /usr/bin/env python

import numpy as n
import pylab as p
import optparse,sys

o  = optparse.OptionParser()

opts,args = o.parse_args(sys.argv[1:])

descriptor = dict(names=('time', 'outT', 'balunT'), formats=('float', 'float', 'float'))
                 
l = 0
for file in args:
    if l == 0:
        data = n.loadtxt(file, dtype=descriptor)
        l += 1
    else:
        data = n.concatenate((data,n.loadtxt(file, dtype=descriptor)))
    
p.plot(data['time'], data['outT'])
p.plot(data['time'], data['balunT'])
p.show()
    

