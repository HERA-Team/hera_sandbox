#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('--k3pk', action='store_true',
    help='Plot Delta^2 instead of P(k)')
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

fq = .16
aa = a.cal.get_aa(opts.cal, n.array([fq]))

uDat, uWgt = {}, {}
for npzfile in args:
    print 'Reading...', npzfile
    dat = n.load(npzfile)
    kpl = dat['kpl']
    keys = dat.files[:]
    keys.remove('kpl')
    
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
    
if opts.k3pk:
    for umag in keys:
        p.loglog(n.abs(kpl),n.abs(kpl)**3*(2*n.pi**2)**-1*n.abs(n.real(uDat[umag])),label=str(umag))
else:
    for umag in keys:
        p.loglog(n.abs(kpl),n.abs(n.real(uDat[umag])),label=str(umag))

p.legend()
p.show()

