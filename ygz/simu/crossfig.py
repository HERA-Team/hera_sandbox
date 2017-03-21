#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
from scipy import signal

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,38)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.03300,'26_46': -0.034, '26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}

dT = 0.033
T0 = 2455700.4
aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
aa.set_jultime(T0)

h = a.healpix.HealpixMap(nside=256)
#h = a.healpix.HealpixMap(nside=64)
#h.set_interpol(False)
#ex,ey,ez = h.px2crd(n.arange(h.map.size), ncrd=3)

bl1x, bl1y, bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
bl2x, bl2y, bl2z = aa.get_baseline(bl2[0],bl2[1],'z')

#fng1=1v

tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
bm = aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/n.abs(tz)   #tz is the Jacobian

bm = n.where(tz > 0.001, bm, 0)
#bm /= bm.sum()
h.map = bm
# C.plot.plot_hmap_ortho(h,mx=0,drng=2)
# p.colorbar()
# p.show()


# #Create equatorial coordinates of the first frame T0
top = n.array([tx,ty,tz], dtype=tx.dtype)

m = n.linalg.inv(aa.eq2top_m)
ex,ey,ez = n.dot(m, top)
eq = n.array([ex,ey,ez], dtype=ex.dtype)

T1 = T0+dT
aa.set_jultime(T1)
m = aa.eq2top_m
t2x,t2y,t2z = n.dot(m, eq)

bl2_prj = t2x*bl2x + t2y*bl2y + t2z*bl2z
bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
fng1 = n.exp(-2j*n.pi*bl1_prj*fq)
fng2 = n.exp(-2j*n.pi*bl2_prj*fq)
# #fng2=1
bm2 = aa[0].bm_response((t2x,t2y,t2z),pol='I')[0]**2#/n.abs(tz)#*n.abs(tzsave)

# #bm = n.ones_like(tx)
# #bm = n.where(tz > 0, bm, 0)
# bm2 = n.where(tz > 0.001, bm2, 0)
# #import IPython; IPython.embed()
# #print bm.sum()
# bm2 /= bm2.sum()
bm_fng2 = bm2 * fng2
bm_fng1 = bm * fng1
h.map = bm_fng1
p.figure()
p.subplot(131)
C.plot.plot_hmap_ortho(h,mode="real")
p.colorbar()
h.map = bm_fng2
p.subplot(132)
C.plot.plot_hmap_ortho(h,mode="real")
p.colorbar()
h.map = bm_fng2*bm_fng1.conj()
p.subplot(133)
C.plot.plot_hmap_ortho(h,mode="real")
p.colorbar()
p.show()

# #IM = a.img.Img(size=200)
# p.ion()
# img = a.img.Img(200,.5)
# tx2,ty2,tz2 = img.get_top((200,200))
# SH = tx2.shape
# tx2,ty2,tz2 = tx2.flatten(),ty2.flatten(),tz2.flatten()
# top2 = n.array([tx2,ty2,tz2], dtype=tx2.dtype)
# ex2,ey2,ez2 = n.dot(m, top2)
# h.map = bm_fng1
# bmfng_proj = h[ex2,ey2,ez2]
# bmfng_proj.shape = SH
# plt = p.imshow(bmfng_proj.imag, vmax=2, vmin=-2,origin='lower')
# p.show()
# #import IPython; IPython.embed()