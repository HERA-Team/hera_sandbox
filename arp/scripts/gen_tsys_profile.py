#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import ephem, sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-m', '--map', dest='map',
    help='Haslam map.')
o.add_option('-c', '--chan', dest='chan', type='int', default=600,
    help='Channel')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
fq = uv['sfreq'] + opts.chan * uv['sdf']
del(uv)

aa = a.cal.get_aa(opts.cal, .001, fq, 1)
im = a.img.Img(size=100, res=.5)
h = a.map.Map(fromfits=opts.map)

tx,ty,tz = im.get_top()
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()

resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[0,0])

lsts_fit = []
dat = []
for ha in n.arange(0,2*n.pi, .1):
    lsts_fit.append(ha)
    print ha
    ex, ey, ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten() 
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = h[ex,ey,ez] * (fq/.408)**-2.52
    tsky.shape = resp.shape
    tsky = n.where(invalid, 0, tsky)
    dat.append(n.sum(tsky * resp) / n.sum(resp))

auto = []
lsts = []
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    uv.select('antennae', 0, 0)
    uv.select('decimate', 2, 0)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        auto.append(d[opts.chan]/aa.passband(i,j))
        aa.set_jultime(t)
        lsts.append(aa.sidereal_time())
    
dat = n.array(dat)
auto = n.array(auto).real
lsts = n.array(lsts)
lsts_fit = n.array(lsts_fit + lsts_fit + lsts_fit)
dat_fit = n.concatenate([dat, dat, dat])
sync_poly = n.polyfit(lsts_fit, dat_fit, deg=12)
sync_auto = n.polyval(sync_poly, lsts)

if True:
    poly = n.polyfit(sync_auto, auto, deg=1)
    print poly
    gain, T_rx = poly
    T_rx /= gain
    print gain, T_rx
else:
    #gain, T_rx = 1057., 92.7
    gain, T_rx = 1143., 76.7

#p.plot(n.arange(0,2*n.pi, .1) * 12 / n.pi, dat)
p.subplot(211)
p.plot(lsts * 12 / n.pi, sync_auto)
p.plot(lsts * 12 / n.pi, auto/gain - T_rx, '.')
p.subplot(212)
p.plot(sync_auto, auto/gain, '.')
Ts = n.arange(0,500,1)
p.plot(Ts, (Ts+T_rx), ':')
p.show()

#resp = a.img.recenter(resp, (100,100))
#tsky = a.img.recenter(tsky, (100,100))

#print n.sum(tsky * resp) / n.sum(resp)

'''
p.subplot(131)
p.imshow(resp)
p.colorbar(shrink=.5)

p.subplot(132)
p.imshow(tsky)
p.colorbar(shrink=.5)

p.subplot(133)
p.imshow(tsky * resp)
p.colorbar(shrink=.5)
'''

p.show()
