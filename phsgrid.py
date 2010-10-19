#!/usr/bin/env python
#
#  phsgrid.py
#  
#
#  Created by Danny Jacobs on 6/17/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
import logging, warnings,os
from mpi4py import MPI
from pylab import *
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,chan=True,pol=True,src=True)
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")
#o.add_option('--N',type='int',
#    help='Number of records to chunk')
opts, args = o.parse_args(sys.argv[1:])



#for each baseline phs and accumulate
#then put the central bit into a matrix


srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]


#Begin: Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#End: Setup MPI

files = comm.scatter(map(list,n.array_split(args,size)))
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])

plot_x = n.zeros([uv['nants'],uv['nants']])
del(uv)
for uvfile in files:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uv, 'cross', opts.pol)
    curtime = 0
    for (uvw,t,(i,j)),d in uv.all():
        if t!= curtime:
            aa.set_jultime(t)
            src.compute(aa)
            curtime = t
        bl = '%d,%d' % (i,j)
        d = aa.phs2src(d, src, i, j)
        flags = n.logical_not(d.mask).astype(n.float)
        gain = n.sqrt(n.average(flags**2))
        ker = n.fft.ifft(flags)        
        d = d.filled(0)
        d = n.fft.ifft(d)
        d, info = a.deconv.clean(d, ker, tol=1e-3)
        d += info['res'] / gain
        d = n.ma.array(d)
        d = n.where(n.isnan(d),0,d)
#        d = n.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], 
#            axis=0)        
        plot_x[i,j] += (n.ma.average(d[:10])+n.ma.average(d[-10:]))/2
del(uv)
plot_xs = comm.gather(plot_x)
if rank!=0: sys.exit()
plot_x = n.zeros_like(plot_xs[0])
for X in plot_xs:
    plot_x += X
n.save('phsgrid',plot_x)
#figure()
#print "plotting"
#matshow(plot_x)
#show()
        