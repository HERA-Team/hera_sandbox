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
a.scripting.add_standard_options(o, ant=True,cal=True,chan=True,pol=True,src=True)
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")
#o.add_option('--N',type='int',
#    help='Number of records to chunk')
opts, args = o.parse_args(sys.argv[1:])

#Logging bits
if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)

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

plot_x = n.zeros([len(aa),len(aa)])
del(uv)
for uvfile in files:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    curtime = 0
    for (uvw,t,(i,j)),d in uv.all():
        if t!= curtime:
            aa.set_jultime(t)
            src.compute(aa)
            curtime = t
        print t,
        sys.stdout.flush()
        flags = n.logical_not(d.mask).astype(n.float)
        gain = n.sqrt(n.average(flags**2))
        ker = n.fft.ifft(flags) 
        if (j-i)>=18: 
            d = n.conjugate(d); 
            print i,j,i-j,'*'
#                    if j%2 != i%2:
#                        d = n.conjugate(d)
#                        print i,j,'*'
        if i%2==1 and j%2==0:
            d = n.conjugate(d); 
            print i,j,i-j,'**'  
#        if (i,j)==(1,30): d = n.conjugate(d); print i,j,'**!'
#        if (i,j)==(1,28): d = n.conjugate(d); print i,j,'**!'
        D = n.ma.copy(d)
        for ni in range(len(aa)):
            for nj in range(ni,len(aa)):                        
                if ni==nj: continue
#                bl = '%d,%d' % (ni,nj)
#                print i,j,ni,nj
#                print n.average(n.abs(D)),
                d = aa.phs2src(D, src, ni, nj)       
#                d = n.where(f, 0, d)
                d = d.filled(0)
                d = n.fft.ifft(d)
                d, info = a.deconv.clean(d, ker, tol=1e-3)
                d += info['res'] / gain
                d = n.ma.array(d)
                d = n.ma.array(n.where(n.isnan(d),0,d))  
#                print n.average(n.abs(d))
        #        d = n.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], 
        #            axis=0)        
                plot_x[ni,nj] += (n.ma.average(d[:5])+n.ma.average(d[-5:]))/2
#                if (ni,nj)==(0,30): print n.average(n.abs(d)),n.max(D),plot_x[ni,nj]
#                sys.stdout.flush()
                del(d)
        plot_xs = comm.gather(plot_x)
        if rank==0:
            plot_x = n.zeros_like(plot_xs[0])
            for X in plot_xs:
                plot_x += X
            n.save('phsgrid',plot_x)
del(uv)
print "%s exiting"%rank

#figure()
#print "plotting"
#matshow(plot_x)
#show()
        