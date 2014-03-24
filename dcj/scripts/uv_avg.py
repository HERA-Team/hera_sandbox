#! /usr/bin/env python
"""
Averages a file into a dict indexed by baseline (i_j)
"""


import aipy as a, numpy as n, sys, os, optparse, pickle
from smooth import smooth

o = optparse.OptionParser()
o.set_usage('uv_avg.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--apply_masknpz',action='store_true',
    help='Applies a mask generated by xrfi_simple.py --npz. Currently applies all mask channels in file named <filename>R.npz')
opts,args = o.parse_args(sys.argv[1:])

for file in args:
    outfile = file+'.avg.pkl'
    print file, ' > ',outfile
    if os.path.exists(outfile):
        print "...exists"
        continue
    uv = a.miriad.UV(file)
    freqs = n.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf'])
    SUM,COUNT = {},{}
    if opts.apply_masknpz and os.path.exists(file+'R.npz'):
        masks = n.load(file+'R.npz')
        M = n.sum(n.array([masks[m] for m in masks.files]),axis=0)
        mask = n.clip(M,0,1)        
    ti = -1
    curtime = 0
    for (uvw,t,(i,j)),d in uv.all():
        if t!= curtime: ti +=1;curtime=t;
        if opts.apply_masknpz: d.mask |= mask[ti,:]
        bl = "%d_%d"%(i,j)
        SUM[bl] =  SUM.get(bl,0) + d.filled(0)
        COUNT[bl] = COUNT.get(bl,0) + n.logical_not(d.mask).astype(int)
    del(uv)
    for k in SUM:
        N = COUNT[k]
        SUM[k][N>0] /= N[N>0]
    SUM['freqs'] = freqs
    SUM['counts'] = COUNT
    F = open(outfile,'w')
    pickle.dump(SUM,F)
    F.close()
    

