#! /usr/bin/env python

import numpy as np
import capo, sys, optparse, aipy

def cov(d1, w1, d2=None, w2=None):
    if d2 is None: d2,w2 = d1.conj(),w1
    d1sum,d1wgt = d1.sum(axis=1), w1.sum(axis=1)
    d2sum,d2wgt = d2.sum(axis=1), w2.sum(axis=1)
    x1,x2 = d1sum / np.where(d1wgt > 0,d1wgt,1), d2sum / np.where(d2wgt > 0,d2wgt,1)
    x1.shape = (-1,1)
    x2.shape = (-1,1)
    d1x = d1 - x1
    d2x = d2 - x2
    C = np.dot(d1x,d2x.T)
    W = np.dot(w1,w2.T)
    return C / np.where(W > 1, W-1, 1)

o = optparse.OptionParser()
aipy.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

info,data,flag = capo.miriad.read_files(args, opts.ant, opts.pol, verbose=True)
covs = {}
for bl in data:
    for pol in data[bl]:
        d,w = data[bl][pol].T, np.logical_not(flag[bl][pol]).astype(np.int).T
        k = str(bl+(pol,))
        covs[k] = cov(d,w)

np.savez('cov.npz', **covs)
#import IPython; IPython.embed()
