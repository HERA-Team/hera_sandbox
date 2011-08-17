#! /usr/bin/env python
"""
Averages a file into a dict indexed by baseline (i_j)
"""


import aipy as a, numpy as n, sys, os, optparse, pickle
from smooth import smooth

o = optparse.OptionParser()
o.set_usage('uv_avg.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for file in args:
    outfile = file+'.avg.pkl'
    print file, ' > ',outfile
    uv = a.miriad.UV(file)
    freqs = n.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf'])
    SUM,COUNT = {},{}
    for (uvw,t,(i,j)),d in uv.all():
        bl = "%d_%d"%(i,j)
        SUM[bl] =  SUM.get(bl,0) + d.filled(0)
        COUNT[bl] = COUNT.get(bl,0) + n.logical_not(d.mask).astype(int)
    del(uv)
    for k in SUM:
        N = COUNT[k]
        SUM[k][N>0] /= N[N>0]
    SUM['freqs'] = freqs
    F = open(outfile,'w')
    pickle.dump(SUM,F)
    F.close()
    

