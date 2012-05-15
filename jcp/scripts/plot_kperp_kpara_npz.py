#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C, optparse, sys

o = optparse.OptionParser()
o.add_option('--fold',action='store_true',
    help='Average negative and positive ks')
opts,args = o.parse_args(sys.argv[1:])

for f in args:
    npz = n.load(f)
    bls = list(npz.files)
    bls.remove('k')
    bls.sort()
    nks = len(npz[bls[0]])
    if opts.fold: nks /= 2
    data = n.zeros((nks,len(bls)))
    p.subplot(211)
    for ind,bl in enumerate(bls):
        if opts.fold: data[:,ind] = (npz[bl][1:nks] + npz[bl][nks:])/2 
        else: data[:,ind] = npz[bl]
        #print npz[bl]
        print f
        p.semilogy(npz['k'],n.abs(npz[bl])/(n.abs(npz['k']**3)),label=bl)
        #p.plot(npz['k'],npz[bl],label=bl)
    p.legend()
    p.subplot(212)
    p.imshow(n.log10(n.abs(data)),aspect='auto',interpolation='nearest',origin='lower')
    p.colorbar()
    p.show()
