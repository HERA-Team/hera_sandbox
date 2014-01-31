#! /usr/bin/env python

import aipy as a, numpy as n, os, sys, glob
import matplotlib.pyplot as p

import optparse
o = optparse.OptionParser()
o.add_option('-l', '--lst', dest='lst', type='float', default=1.0,
        help='Which lst bin to use for fractional spread')
a.scripting.add_standard_options(o,chan=True)

opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan,uv['nchan'])

width_data = {}
lst = []

for file in args:
    print file
    curtime = None
    uv = a.miriad.UV(file)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            curtime = t
            lst = uv['lst']
            d = d.real
            if n.abs(lst - opts.lst) < 0.01:
                print lst
                if not width_data.has_key(int(t)): 
                    width_data[int(t)] = [d.take(chans)] 
                else: 
                    width_data[int(t)].append(d.take(chans))

full_data = []
for keys in width_data.keys():
    full_data.append(n.sum(width_data[keys])/len(width_data[keys]))
full_data = n.array(full_data)    
index = n.where(full_data != 0.0)
mean = n.mean(full_data[index])
max = n.max(full_data[index])
min = n.min(full_data[index])
p.plot(full_data[index])
print max,min,mean

print 'fractional spread = %f'%((max - min)/mean)
                
p.show()
            
                

