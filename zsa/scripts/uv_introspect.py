#!/usr/bin/env python
import aipy as a
import sys

for file in sys.argv[1:]:
    uv = a.miriad.UV(file)
    ants = []
    times = []
    lst  = []
    curtime = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            times.append(t)
            lst.append(uv['lst'])
            curtime = t 
        if '%d_%d'%(i,j) in ants:continue
        else: ants.append('%d,%d'%(i,j))
    
    print 'Summary for %s : '%file 
    print 50*'-'
    print '\t Antennas ( %d baselines) : '%len(ants), ants
    print '\t Time range %f to %f (%d samples in %f increments)'%(times[0],times[-1],len(times), times[2]-times[1])
    print '\t LST range %f to %f (%dsamples in %f increments)'%(lst[0],lst[-1],len(lst), lst[2]-lst[1])
    print '\n\n'


    
    
         
