#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('--ubin', action='store_true',
    help='Bin baselines by u-magnitude')
o.add_option('--k3pk', action='store_true',
    help='Plot Delta^2 instead of P(k)')
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

npzfile = args[0]
dat = n.load(npzfile)
keys = dat.files[:]
keys.remove('kpl')

fq = .16
aa = a.cal.get_aa(opts.cal, n.array([fq]))


if opts.ubin:
    uDat, uWgt, blDat = {}, {}, {}
    #bin in umag
    for bl in keys:
        if bl[0] == 'w': continue
        wbl = 'w'+str(bl)
        i,j = a.miriad.bl2ij(bl)
        crd = aa.get_baseline(i,j)*fq
        umag = (crd[0]**2 + crd[1]**2)**.5
        umag = str(2**int(n.around(n.log2(umag.clip(0.5,n.Inf)))))
        uDat[umag] = uDat.get(umag,0) + dat[bl]
        uWgt[umag] = uWgt.get(umag,0) + dat[wbl]

    keys = uDat.keys()
    for umag in keys:
        uDat[umag] /= uWgt[umag]
else:
    uDat = dat

if opts.k3pk:
    for umag in keys:
        p.loglog(n.abs(dat['kpl']),n.abs(dat['kpl'])**3*(n.pi*2)**-2*n.abs(n.real(uDat[umag])),label=str(umag))
else:
    for umag in keys:
        p.loglog(n.abs(dat['kpl']),n.abs(n.real(uDat[umag])),label=str(umag))

p.legend()
p.show()

