#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,38)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.032557,'26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}

try: ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
except(KeyError): ver = 0.
print 'DelT = ', ver
TT = n.arange(2455700,2455701,1/a.const.sidereal_day*42.9*5) #*5 for speed
T0 = 2455700.5
T1 = T0+ver
#############################################################################
aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
#h = a.healpix.HealpixMap(nside=128)
h = a.healpix.HealpixMap(nside=64)
tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
bl1x, bl1y, bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
bl2x, bl2y, bl2z = aa.get_baseline(bl2[0],bl2[1],'z')

bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
fng1 = n.exp(-2j*n.pi*bl1_prj*fq)
bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2
bm = n.where(tz > 0, bm, 0)
print bm.sum()
bm /= bm.sum()
bm_fng1 = bm * fng1


top = n.array([tx,ty,tz], dtype=tx.dtype)

aa.set_jultime(T0)
m = n.linalg.inv(aa.eq2top_m)
ex,ey,ez = n.dot(m, top)
eq = n.array([ex,ey,ez], dtype=ex.dtype)

aa.set_jultime(T1)
m = aa.eq2top_m
tx,ty,tz = n.dot(m, eq)

bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
fng2 = n.exp(-2j*n.pi*bl2_prj*fq)
bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2
bm = n.where(tz > 0, bm, 0)
print bm.sum()
bm /= bm.sum()
bm_fng2 = bm * fng2

#import IPython; IPython.embed()
print bm_fng1.shape, bm_fng2.shape
print n.sum(bm_fng1.conj()*bm_fng2)
import IPython; IPython.embed()