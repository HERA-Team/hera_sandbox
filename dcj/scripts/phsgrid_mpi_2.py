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
o.add_option('--exp',action='store_true',
    help='Munge the data experimentally.. kinda thing')
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


def phsdlyavg(data):
    comm.bcast(2,root=0)
    data = comm.scatter([0]*1+n.array_split(data,size-1),root=0)
    plot_xs = comm.gather([0]*1)
    return n.sum(plot_xs[1:],axis=0)

#for each baseline phs and accumulate
#then put the central bit into a matrix


srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]

rewire_run1 = { 
    0: {'x':0, 'y':1},
    1: {'x':2, 'y':3},
    2: {'x':4, 'y':5},
    3: {'x':6, 'y':7},
    4: {'x':8, 'y':9},
    5: {'x':10, 'y':11},
    6: {'x':12, 'y':13},
    7: {'x':14, 'y':15},
    8: {'x':16, 'y':17},
    9: {'x':18, 'y':19},
    10: {'x':20, 'y':21},
    11: {'x':22, 'y':23},
    12: {'x':24, 'y':25},
    13: {'x':26, 'y':27},
    14: {'x':28, 'y':29},
    15: {'x':30, 'y':31},
}

#Begin: Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#End: Setup MPI

#files = comm.scatter(map(list,n.array_split(args,size)))
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])

uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                     ('t', n.float), 
                     ('bl', '<i4', 2)]),
                ('spec', n.complex64,uv['nchan']),
                ('mask',n.bool,uv['nchan'])])
plot_x = n.zeros([len(aa),len(aa)])
plot_x0 = n.zeros_like(plot_x)
del(uv)
if rank==0:
    print "making a per baseline phase diagram across: %i nodes"%size
    if opts.exp: print "Performing some dumb experiment"
    for uvfile in args:
        print 'Reading', uvfile
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, 'cross', 'xx,yy,xy,yx')
        curtime = 0
        data = n.array([],dtype=uv_aipy,ndmin=1)
        for (uvw,t,(i,j)),d in uv.all():
            if opts.exp:
                t -= 2.0/24
                pi,pj = a.miriad.pol2str[uv['pol']]
                print pi,pj
                if i==j and (pi,pj)==('y','x'): continue
                if (j-i)>16:
                     d = n.conjugate(d)
                     if pi != pj: 
                         pi,pj = pj,pi
                         d = n.conjugate(d)
                if (pi,pj)==('y','x'): d = n.conjugate(d)
                i = rewire_run1[i][pi]
                j = rewire_run1[j][pj]
                if i> j:
                    i,j=j,i
                    d = n.conjugate(d)
#                if (j-i)>=18: 
#                    d = n.conjugate(d); 
#                    print i,j,i-j,'*'
#                    if j-i == 19 and i%2==0 and j%2==1: pass
#                    elif i%2 != j%2:
#                        print i,j,'-->',                    
#                        i += (-1)**(i%2)
#                        j += (-1)**(j%2)
#                        print i,j,'*'
#                        d = n.conjugate(d)
##                        print i,j,'*'
#                if i%2==1 and j%2==0:
#                    d = n.conjugate(d); 
##                    print i,j,i-j,'**'                    
                sys.stdout.flush()
            p = (uvw,t,(i,j))
            if t!= curtime:
                if t-curtime<2e5:#do the scatter and calc from last time step (if its not the first time step)
                    plot_x0 += phsdlyavg(data)
                    data = n.array([],dtype=uv_aipy,ndmin=1)
                aa.set_jultime(t)
                src.compute(aa)
                curtime = t
            mask = d.mask
            d = d.data
            rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
            data = n.concatenate((data,rec))
        plot_x0 += phsdlyavg(data)
        n.save('phsgrid_mpi_2',plot_x0)

    del(uv)
    comm.bcast(0)
else:
    while(True):
        cmd = comm.bcast(root=0)
        if cmd==2:
            plot_x = n.zeros_like(plot_x)
            data = comm.scatter(root=0)
#            print '.',
#            sys.stdout.flush()
            curtime = 0
            for p,d,f in data:
                i,j = p['bl']
                t = p['t']
                if t!=curtime:
                    aa.set_jultime(t)
                    src.compute(aa)
                    curtime=t
                d = aa.phs2src(d, src, i, j)
                flags = n.logical_not(f).astype(n.float)
                gain = n.sqrt(n.average(flags**2))
                ker = n.fft.ifft(flags)        
                d = n.where(f, 0, d)
                d = n.fft.ifft(d)
                d, info = a.deconv.clean(d, ker, tol=1e-3)
                d += info['res'] / gain
                d = n.ma.array(n.where(n.isnan(d),0,d))
                d = n.ma.abs(n.ma.masked_less_equal(d,0))
                plot_x[i,j] += (n.ma.average(d[:5])+n.ma.average(d[-5:]))/2
                plot_x = n.where(n.isnan(plot_x),0,plot_x)
            comm.gather(plot_x)
            plot_x = n.zeros_like(plot_x)
        elif cmd==0:
            sys.exit()

#figure()
#print "plotting"
#matshow(plot_x)
#show()
        