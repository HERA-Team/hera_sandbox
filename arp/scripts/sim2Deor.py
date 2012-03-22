#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

p.ion()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv)

#aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa = a.cal.get_aa(opts.cal, n.array([.1738]))

h = a.healpix.HealpixMap(nside=64)
#h = a.healpix.HealpixMap(nside=32)
h.map = n.random.normal(size=h.map.size)
h.set_interpol(True)

#im = a.img.Img(size=200, res=.5)
SZ = 200
im = a.img.Img(size=SZ, res=.5)
tx,ty,tz = im.get_top()
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()

resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[0,0])

lsts = []
uvdata0,uvdata1,uvdata2,uvdata3 = [],[],[],[]
PLOT = None
for ha in n.arange(-n.pi/16, n.pi/16, 2*n.pi / (24*60)):
    print ha
    lsts.append(ha)
    ex,ey,ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten()
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    Tsky = h[ex,ey,ez]
    Tsky.shape = resp.shape
    img = Tsky * resp
    img = a.img.recenter(img, (img.shape[0]/2, img.shape[1]/2))
    uv = n.fft.fft2(img)
    plt_uv0 = a.img.recenter(n.abs(uv), n.array(uv.shape)/2)[SZ:SZ+64,SZ:SZ+64]
    plt_uv1 = a.img.recenter(n.angle(uv), n.array(uv.shape)/2)[SZ:SZ+64,SZ:SZ+64]
    if PLOT is None:
        PLOT = {}
        p.subplot(221)
        PLOT['img0'] = p.imshow(img, vmax=img.max(), vmin=-img.max(), origin='lower')
        p.subplot(222)
        PLOT['img1'] = p.imshow(n.abs(img), vmax=img.max(), vmin=0, origin='lower')
        p.subplot(223)
        PLOT['uv0'] = p.imshow(plt_uv0, origin='lower')
        p.subplot(224)
        PLOT['uv1'] = p.imshow(plt_uv1, origin='lower')
    else:
        PLOT['img0'].set_data(img)
        PLOT['img1'].set_data(n.abs(img))
        PLOT['uv0'].set_data(plt_uv0)
        PLOT['uv1'].set_data(plt_uv1)
    uvdata0.append(uv[0,:uv.shape[1]/2])
    uvdata1.append(uv[16,:uv.shape[1]/2])
    uvdata2.append(uv[32,:uv.shape[1]/2])
    uvdata3.append(uv[:uv.shape[0]/2,0])
    p.draw()

uvdata0 = n.array(uvdata0); uvdata0 /= n.abs(uvdata0)
uvdata1 = n.array(uvdata1); uvdata1 /= n.abs(uvdata1)
uvdata2 = n.array(uvdata2); uvdata2 /= n.abs(uvdata2)
uvdata3 = n.array(uvdata3); uvdata3 /= n.abs(uvdata3)
lsts = n.array(lsts)

p.ioff()
p.clf()
nacc = n.arange(uvdata0.shape[0])+1
for step in [8,16,32,64]:
    p.subplot(2,2,1)
    p.loglog(nacc*60, n.abs(n.cumsum(uvdata0[:,step]) / nacc), label='uv=%2d, 0'%(step*.5))
    p.subplot(2,2,2)
    p.loglog(nacc*60, n.abs(n.cumsum(uvdata1[:,step]) / nacc), label='uv=%2d, 8'%(step*.5))
    p.subplot(2,2,3)
    p.loglog(nacc*60, n.abs(n.cumsum(uvdata2[:,step]) / nacc), label='uv=%2d,16'%(step*.5))
    p.subplot(2,2,4)
    p.loglog(nacc*60, n.abs(n.cumsum(uvdata3[:,step]) / nacc), label='uv= 0,%2d'%(step*.5))

for i in range(4):
    p.subplot(2,2,i+1)
    p.loglog(nacc*60, 3*nacc**-.5, 'k:')
    p.xlim(1e2, 1e4)
    p.ylim(1e-2, 1e1)
    p.legend(loc='best')
    p.grid()

p.show()
