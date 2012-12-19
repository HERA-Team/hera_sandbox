#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

aa = a.cal.get_aa('psa898_v002', n.array([.150]))
im = a.img.Img(size=200, res=.5)
tx,ty,tz = im.get_top(center=(200,200))
valid = n.logical_not(tx.mask)
tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
amp = aa[0].bm_response((tx,ty,tz),pol='x')**2
bx,by,bz = aa.get_baseline(0,1,'z')
print [bx,by,bz]
for fq in n.arange(.1,.2,.01):
    phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz))
    phs.shape = amp.shape = im.uv.shape
    amp = n.where(valid, amp, 0)
    phs = n.where(valid, phs, 0)

    dc_response = n.sum(amp*phs)/n.sum(amp)
    print fq, dc_response, n.real(dc_response * 100 / C.pspec.jy2T(fq))
    if False:
        p.subplot(121); p.imshow(amp)
        p.subplot(122); p.imshow(amp*phs.real,vmax=1,vmin=-1)
        p.show()

