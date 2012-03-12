#!/usr/bin/env python
#
#  d_ogram.py
#  
#
#  Created by Danny Jacobs on 10/7/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,ant=True,pol=True,dec=True,chan=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

D = []
for file in args:
    uv = a.miriad.UV(file)
    aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
    aa.select_chans(chans)
    
    for (uvw,t,(i,j)),dat in uv.all():
#        bl = a.miriad.ij2bl(i,j)
        if i==j:continue
        uvw = aa.gen_uvw(i,j).squeeze()
        bl = n.sqrt(n.dot(aa.gen_uvw(i,j).squeeze().T,aa.gen_uvw(i,j).squeeze()))[0,:]
#        print uvw.shape,bl.shape,dat.take(chans).shape
        D.append(zip(bl,dat.take(chans)))

S = n.array(D)[:,:,1]
D = n.concatenate(D)
C = n.argwhere(n.abs(D[:,1])>10).squeeze()
print C
D = D[C]
#print n.max(D[:,0]),n.max(D[:,1]),D
p.subplot(311)
p.semilogy(D[:,0],n.abs(D[:,1]),'.',alpha=0.1,markersize=1)
p.ylabel('raw visibility')
p.xlabel('baseline length [wl]')
p.subplot(312)
print n.sqrt(n.max(D.shape))
print n.max(n.log10(n.abs(D[:,1]))),n.min(n.log10(n.abs(D[:,1])))
p.hist(n.log10(n.abs(D[:,1])),bins=n.sqrt(n.max(D.shape)),histtype='step',log=1)
p.xlabel('raw visibility amplitude')
p.subplot(313)
for d in S:
    p.semilogy(aa.get_afreqs(),n.ma.abs(d),'.',color='b',alpha=0.1,markersize=3)
    p.xlabel('freq [GHz]')
p.show()