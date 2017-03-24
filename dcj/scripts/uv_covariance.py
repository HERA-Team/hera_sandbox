#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p, capo, capo.frf_conv as fringe
import glob, optparse, sys, random
from capo.dcj import condition_goodness
from matplotlib.pyplot import *
def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex128)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1) #normalization
    return (n.dot(X, X.T.conj()) / fact).squeeze()



o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True,dec=True)
o.add_option('--plot',action='store_true',
    help='Plot plots')
opts,args = o.parse_args(sys.argv[1:])

print "loading data..",
info, data, flag = capo.miriad.read_files(args,opts.ant,opts.pol,decimate=opts.decimate,recast_as_array=True)
print "[finished]"
print "found {n} times".format(n=len(info['times']))
Cs = []
conds = []
goods = []
chans = None
for bl in data:
    for pol in data[bl]:
        if chans is None:
            chans = a.scripting.parse_chans(opts.chan,data[bl][pol].shape[1])
        D = data[bl][pol][:,chans].T.astype(np.complex128)
        C = cov(D)
        U,S,V = n.linalg.svd(C.conj())
        _C = n.einsum('ij,j,jk', V.T, 1./S, U.T)
        #show()
        Cs.append(C)
        conds.append(n.log(n.abs(n.linalg.cond(C))/n.log(2)))
        goods.append(condition_goodness(C))
        print pol,bl,n.round(conds[-1],2),goods[-1],D.dtype
        if opts.plot:
            subplot(321)
            imshow(n.log(n.abs(C)))
            subplot(322)
            imshow(n.log(n.abs(_C)))
            subplot(312)
            semilogy(S)
            subplot(313)
            imshow(n.abs(D),interpolation='nearest',aspect='auto')
            tight_layout()
            show()
            
            


