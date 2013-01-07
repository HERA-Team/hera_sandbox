#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('--k3pk', action='store_true',
    help='Plot Delta^2 instead of P(k)')
o.add_option('--umax', type='float', default=n.Inf,
    help='Only show baselines shorter than this value')   
o.add_option('--umin', type='float', default=0.,
    help='Only show baselines longer than this value')   
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
        if i == 40 or j == 40: continue
        if i == 55 or j == 55: continue
        crd = aa.get_baseline(i,j)*fq
        umag = (crd[0]**2 + crd[1]**2)**.5
        #pick predominantly east-west baselines
        #if crd[0]**2 < .85 * umag**2: continue
        if umag > opts.umax: continue
        if umag < opts.umin: continue
        umag = str(2**int(n.around(n.log2(umag.clip(0.5,n.Inf)))))
        uDat[umag] = uDat.get(umag,0) + dat[bl]
        uWgt[umag] = uWgt.get(umag,0) + dat[wbl]
    
keys = uDat.keys()
for umag in keys:
    uDat[umag] /= uWgt[umag]
    

nbin = len(keys)
for ind,umag in enumerate(keys):
    p.subplot(n.ceil(n.sqrt(nbin)),n.ceil(n.sqrt(nbin)),ind+1)
    if opts.k3pk: f = n.abs(kpl)**3*(2*n.pi**2)**-1
    else: f = 1.
    noise_ind = n.where(n.abs(kpl > 1.))
    noise = n.mean(n.abs(n.real(uDat[umag][noise_ind])))
    noise *= n.ones_like(kpl)
    p.loglog(n.abs(kpl),f*noise,color='k',lw=3)
    #fg = n.array([C.pspec.dk_du(C.pspec.f2z(fq)) * n.float(umag),C.pspec.dk_du(C.pspec.f2z(fq)) * n.float(umag)])
    fg = n.array([C.pspec.dk_deta(C.pspec.f2z(fq))*n.float(umag)/fq,C.pspec.dk_deta(C.pspec.f2z(fq))*n.float(umag)/fq])
    yfg = n.array([1e-5,1e9])
    label = str(umag)
    p.loglog(fg,yfg,color='r',lw=3)
    p.loglog(n.abs(kpl),f*n.abs(n.real(uDat[umag])),label=label,color='b') 
    p.ylim(1e-1,1e7)
    p.legend()
p.show()

