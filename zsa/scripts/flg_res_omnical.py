#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

seps = [6,11,12,13,21,24,25,26,28,30,32,35,44,47,48,49,50,52,64,84,85,99]
#seps = [24]
seps = ['sep%d'%i for i in seps]

sig = {}
for sep in seps:
    d,w = [],[]
    for filename in sys.argv[1:]:
        print 'Reading', filename, sep
        npz = n.load(filename)
        d.append(npz[sep])
        _wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
        wgt = npz['wgt_'+sep]
        w.append(wgt * _wgt)
    d,w = n.concatenate(d),n.concatenate(w)
    w = n.where(n.isnan(d), 0, w)
    w = n.where(w < 1e-8, 0, w)

    d = n.where(w > 0, d, 0)
    print n.any(n.isnan(d))

    sig[sep] = n.sqrt(n.median(n.abs(d)**2, axis=0))
    sig[sep].shape = (1,sig[sep].size)
    #w = n.where(n.abs(d) > 3*sig, 0, w)
    #d = n.where(w > 0, d, 0)

for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    info = dict(npz)
    _wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
    for sep in seps:
        d = npz[sep]
        w = npz['wgt_'+sep]
        w = n.where(n.isnan(d), 0, w)
        w = n.where(_wgt < 1e-8, 0, w)
        w = n.where(n.abs(d) > 3*sig[sep], 0, w)
        d = n.where(w > 0, d, 0)
        print n.all(w == npz['wgt_'+sep])
        info[sep] = d
        info['wgt_'+sep] = w
    outfile = filename.split('.npz')[-2]+'f.npz'
    print 'Writing', outfile
    n.savez(outfile, **info)
        
