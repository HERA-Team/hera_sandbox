#! /usr/bin/env python
import numpy as n
import sys

dly = {}
gain = {}
times = []
for filename in sys.argv[1:]:
    print 'Reading', filename
    f = open(filename)
    npz = n.load(f)
    C_phs = npz['C_phs']
    C_amp = 10**npz['C_amp']
    antpos = npz['antpos']
    time = npz['time']
#    pols = npz['pols']
    pols = ['xx']
    times.append(time)
    for n,phs in enumerate(C_phs.flatten()):
        if n%8==0:print('\n')
        print '%d7.5'%phs,5*'',
        
#    for pi, pol in enumerate(pols):
#        if not dly.has_key(pi): dly[pi],gain[pi] = {},{}
#        for i,tau,g in zip(antpos.flatten(), C_phs[pi].flatten(), C_amp[pi].flatten()):
#            dly[pi][i] = dly[pi].get(i,[]) + [tau]
#            gain[pi][i] = gain[pi].get(i,[]) + [g]
    f.close()
