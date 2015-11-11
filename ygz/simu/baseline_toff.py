#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
#TSKY,TGND = 300,100 # K
TSKY,TGND = 1,0 # K

aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
#h = a.healpix.HealpixMap(nside=128)
h = a.healpix.HealpixMap(nside=64)
tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
bl1x, bl1y, bl1z = aa.get_baseline(0,26,'z')
bl2x, bl2y, bl2z = aa.get_baseline(0,38,'z')
bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
fng1 = n.exp(2j*n.pi*bl1_prj*fq)
fng2 = n.exp(2j*n.pi*bl2_prj*fq)
bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2
bm = n.where(tz > 0, bm, 0); bm /= bm.sum()
bm_fng1 = bm * fng1
bm_fng2 = bm * fng2

#if True:
#    img = a.img.Img(200,.5)
#    tx,ty,tz = img.get_top((200,200))
#    SH = tx.shape
top = n.array([tx,ty,tz], dtype=tx.dtype)
#top.shape = (3,-1)
#import IPython; IPython.embed()
#for lst in n.arange(0,2*n.pi,2*n.pi/a.const.sidereal_day*42.9):
#plt = None
corr = []
TT = n.arange(2455700.5,2455701.5,1/a.const.sidereal_day*42.9*5)
for i in xrange(3):
    print i
    sky = n.random.normal(size=h.map.size)
    h.map = sky # assume sky is in eq coord
    vis1,vis2 = [],[]
    aa.set_jultime(2455700)
    m = n.linalg.inv(aa.eq2top_m)
    ex,ey,ez = n.dot(m, top)
    sky_prj1 = h[ex,ey,ez]
    for jd in TT:
        # convert tx,ty,tz to ex,ey,ez (time dependent)
        #print jd
        aa.set_jultime(jd)
        m = n.linalg.inv(aa.eq2top_m)
        ex,ey,ez = n.dot(m, top)
        sky_prj2 = h[ex,ey,ez]
        #sky_prj.shape = SH
        #if plt is None: plt = p.imshow(sky_prj, vmax=2, vmin=-2)
        #else:
        #    plt.set_data(sky_prj)
        #    p.draw()
        vis1.append(n.sum(sky_prj1 * bm_fng1))
        vis2.append(n.sum(sky_prj2 * bm_fng2))
        
    vis1,vis2 = n.array(vis1), n.array(vis2)

    _vis1,_vis2 = n.fft.ifft(vis1), n.fft.ifft(vis2)
    corr.append(n.fft.fft(_vis1*n.conj(_vis2)))
    #p.plot(n.abs(corr[i]))
corr = n.array(corr)
print corr.shape
#p.show()
p.plot(TT,n.abs(n.mean(corr,axis=0)))
import IPython; IPython.embed()
