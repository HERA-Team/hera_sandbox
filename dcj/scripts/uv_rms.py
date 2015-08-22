#! /usr/bin/env python
"""
Takes derivatives of input uv files and spits out an rms spectrum (npz file).

usage: uv_rms.py -f -t -p *uv

Limitations: only works on one pol.

"""

import aipy as a, numpy as n, sys, optparse, os 
from capo.pspec import jy2T
from pylab import *
from capo.arp import get_dict_of_uv_data
o = optparse.OptionParser()
o.set_usage('uv_rms.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True,cal=True,pol=True)
o.add_option('-f', action='store_true',default=False,
    help='Do freq derivative.')
o.add_option('-t',action='store_true',default=False,
    help='Do time derivative.')

o.add_option('-i',action='store_true',
    help='clear saved npz and rerun on data')
o.add_option('-v',action='store_true',
    help='Verbose. print more debug stuff')
o.add_option('--bl',action='store_true',
    help='Difference across baselines!.')
o.add_option('--label',type=str,default='None',
    help='add a label to the output data file. useful as a legend later')
o.add_option('--bl_avg',action='store_true',
    help='First, average across baselines!.')

    
opts, args = o.parse_args(sys.argv[1:])
if opts.cal != None:
    uv = a.miriad.UV(args[0])
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    del(uv)
else: aa = None
uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan']).squeeze()
print len(freqs)
if opts.f:
    freqs = freqs[:-1]
del(uv)
def intsquare(x):
    return int(n.sqrt(x)),n.ceil(x/int(n.sqrt(x)))

if len(opts.pol.split(','))>1:
    print "error invalid pol input"
    print "please choose one pol to rms"
    print "exiting"
    sys.exit(1)
if not os.path.exists('uv_rms.npz') or opts.i:
    t,dat,flg = get_dict_of_uv_data(args,opts.ant,opts.pol,verbose=opts.v)
    bls = dat.keys()
    pol = opts.pol
    D = n.ma.masked_where([flg[bl][pol] for bl in bls],[dat[bl][pol] for bl in bls])
    print D.shape
    W = 1 - n.array([flg[bl][pol].astype(n.float) for bl in bls])
    if opts.bl_avg:
        D = n.ma.sum(D,axis=0)
        W = n.ma.sum(W,axis=0)
        D[W>0] /= W[W>0]
        D.shape = (1,) + D.shape
        W.shape = (1,)+W.shape
    #NB: differencing picks up an extra factor of root 2, divide it out for each derivative
    #see Eq 4.16 of Jacobs thesis
    if opts.bl and not (D.shape[0]%2):
        D = n.ma.diff(D,axis=0)/n.sqrt(2)
        if opts.v:
            print "diffing baselines"
    if opts.t:
        D = n.ma.diff(D,axis=1)/n.sqrt(2)
    if opts.f:
        D = n.ma.diff(D,axis=2)/n.sqrt(2)
        W = W[:,:,:-1]#if we're diffing by freqs drop the topmost channel flag

    var = n.abs(n.ma.sum(n.ma.sum(D*n.conj(D),axis=0),axis=0)) #bl=0,time=1
    w = n.sum(n.sum(W,axis=0),axis=0) #bl=0,time=1
    Trms         = var
    Trms[w>0]   /= w[w>0]
    Trms         = n.sqrt(Trms)*jy2T(freqs) 

    if opts.v:
        mypoint = w.argmax()
        print var[mypoint]
        print w[mypoint]
        print Trms[mypoint]

    print "writing uv_rms.npz"
    n.savez('uv_rms.npz',freqs=freqs,weights=w,Trms=Trms.filled(0),mask=Trms.mask,label=opts.label)
    if opts.v:
        rows,cols = intsquare(D.shape[0])
        for i in xrange(D.shape[0]):
            subplot(rows,cols,i+1)
            title(a.miriad.bl2ij(bls[i]))
            imshow(n.abs(D[i])*jy2T(freqs),aspect='auto',vmin=0,vmax=30)
            colorbar()
        figure()
        subplot(211)
        V = n.sqrt(n.abs(n.ma.sum(D*n.conj(D),axis=0)/n.sum(W,axis=0)*jy2T(freqs)))
        imshow(V,aspect='auto',vmin=0,vmax=30)
        colorbar()
        subplot(212)
        semilogy(freqs,Trms)
        show()
