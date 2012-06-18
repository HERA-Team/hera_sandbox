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
    
    #bin in uv-plane
    sums = []
    for bl in keys:
        sums.append(n.sum(dat[bl]))
    for bl in keys:
        if bl[0] == 'w': continue
        wbl = 'w'+str(bl)
        i,j = a.miriad.bl2ij(bl)
        if i == 40 or j == 40: continue
        if i == 55 or j == 55: continue
        crd = aa.get_baseline(i,j)*fq
        u,v = crd[0],crd[1]
        bin = C.pspec.uv2bin(u,v,0)
        #try flagging on sum bin per 10 min file
        #if n.sum(dat[bl]) > 3*n.std(sums): flag = 0
        #else: flag = 1
        flag = 1
        uDat[bin] = uDat.get(bin,0) + flag * dat[bl]
        uWgt[bin] = uWgt.get(bin,0) + flag * dat[wbl]

slices = ['sum','hor-min','hor','hor-plus']

size = 75
center = size/2
uvplane = {}
uvplane_wgt = {} #n.zeros((size,size))
keys = uDat.keys()
for slice in slices:
    print slice
    uvplane[slice] = n.zeros((size,size))
    uvplane_wgt[slice] = n.zeros((size,size))
    vcnt = 0.
    for bin in keys:
        u,v,lst = C.pspec.bin2uv(bin)
        umag = n.sqrt(u**2 + v**2)
        iu,iv = n.round(.25*u+center),n.round(.25*v+center)
        #hor = C.pspec.dk_deta(C.pspec.f2z(fq))*n.float(umag)/fq
        hor = .15
        kpl = n.abs(kpl)
        if slice == 'sum':
            uDat[bin] /= uWgt[bin]
            valid = n.ones_like(kpl)
        if slice == 'hor':
            lo,hi = .75*hor,1.5*hor
            valid = n.where(kpl > lo,1,0)*n.where(kpl < hi,1,0)
        if slice == 'hor-min':
            lo,hi = 0*hor,.75*hor
            valid = n.where(kpl > lo,1,0)*n.where(kpl < hi,1,0)
        if slice == 'hor-plus':
            lo,hi = 1.5*hor,3*hor
            valid = n.where(kpl > lo,1,0)*n.where(kpl < hi,1,0)
        vcnt += n.sum(valid)
        uvplane[slice][iu,iv] += n.sum(uDat[bin].compress(valid))
        uvplane_wgt[slice][iu,iv] += 1.
    print vcnt

xplots = n.ceil(n.sqrt(len(slices)))
yplots = n.floor(n.sqrt(len(slices)))
vmax = [7.,5.,5.,5.]
vmin = [4.,0.,0.,0.]
for plot,slice in enumerate(slices):
    uvplane[slice] /= uvplane_wgt[slice].clip(1.,n.Inf)
    p.subplot(xplots,yplots,plot+1)
    p.title(slice)
    #p.imshow(n.log10(uvplane[slice].clip(1.,n.Inf)),interpolation='nearest',vmax=vmax[plot],vmin=vmin[plot])
    p.imshow(n.log10(n.abs(uvplane[slice])),interpolation='nearest',vmax=vmax[plot],vmin=vmin[plot])
    #p.imshow(uvplane_wgt[slice].clip(0.,1.),interpolation='nearest')
    #p.colorbar()
p.show()
