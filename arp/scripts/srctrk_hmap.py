#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, re
import glob, os
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_bm_hmap.py [options] *.hmap')
opts,args = o.parse_args(sys.argv[1:])

srcfiles = sys.argv[1:]

def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / dat.size
def beamsq_area(dat):
    return n.sum(dat**2) * 4 * n.pi / dat.size

def fng_wgt(fng_bin, fq, bl_ew_len, cen164=.000915):
    # This is filter used in fringe_rate_filter_weighted
    cen = cen164 * (fq * bl_ew_len / 16.43) # Calibrated to .164 GHz, 30m baseline
    wid = .000160 * (fq * bl_ew_len / 16.43) # Calibrated to .164 GHz, 30m baseline
    return n.exp(-(fng_bin-cen)**2/(2*wid**2))

RES = .01
SZ = (int(n.pi/2/RES), int(2*n.pi/RES)+2)
alt,az = n.indices(SZ)
alt = n.pi/2 - alt.astype(n.float) * RES
az = az.astype(n.float) * RES
x,y,z = a.coord.azalt2top((az.flatten(), alt.flatten()))
m = Basemap(projection='ortho', lat_0=90, lon_0=180)
cx,cy = m(az * a.img.rad2deg, alt * a.img.rad2deg)

def plot_hmap(hmaps, max=0, min=-2):
    data = hmap[x,y,z]
    data = n.log10(data).clip(min,max)
    data.shape = az.shape
    m.drawmapboundary()
    m.drawmeridians(n.arange(0,360,30))
    m.drawparallels(n.arange(0,90,10))
    step = (max - min) / 35.
    levels = n.arange(min-step, max+step, step)
    m.contourf(cx, cy, data, levels)
    #p.colorbar(shrink=.5)
    return m

hmaps = [a.map.Map(nside=nside) for nside in [64,32,16,8]]

for filename in srcfiles:
    print '    Reading', filename
    npz = n.load(filename)
    fq = npz['freq']
    print '    Freq:', fq
    sdat = n.abs(npz['spec'])
    valid = n.where(sdat == 0, 0, 1)
    sdat = sdat.compress(valid)
    sx = npz['x'].compress(valid)
    sy = npz['y'].compress(valid)
    sz = npz['z'].compress(valid)
    swgt = n.ones_like(sdat)
    for hmap in hmaps: hmap.add((sx,sy,sz), swgt, sdat)

for hmap in hmaps: hmap.reset_wgt()
d,w = n.zeros(hmaps[0].npix()), n.zeros(hmaps[0].npix())
mx,my,mz = hmaps[0].px2crd(n.arange(hmaps[0].npix()))
for hmap in hmaps:
    wgt = hmap.nside()
    w += hmap.wgt[mx,my,mz] * wgt
    d += hmap.map[mx,my,mz] * wgt
data  = n.where(w > 0, d/w, 0)
hmap = hmaps[0]
hmap.map.map,hmap.wgt.map = d,w
#p.subplot(121)
plot_hmap(hmap)
p.colorbar(shrink=.5)

Omega_P, Omega_PP = beam_area(data), beamsq_area(data)
print 'Omega_P:', Omega_P
print 'Omega_PP:', Omega_PP
print 'Omega_eff:', Omega_P**2/Omega_PP
bl_ew_len = 100. # ns
fng_mx = (fq * bl_ew_len * 2*n.pi / a.const.sidereal_day)
fng_mn = -0.5 * fng_mx
fng_bins = n.linspace(fng_mn, fng_mx, 1000)
fng_wgts = fng_wgt(fng_bins, fq, bl_ew_len)
t_eff = 1. / n.average(fng_wgts**2)
signal_attenuation = Omega_PP / 0.323
print 't_fng/t_sky:', t_eff
print 'signal attenuation:', signal_attenuation
print 'Delta^2_N,fng / Delta^2_N,sky single mode:', 1. / signal_attenuation / t_eff
print 'Delta^2_N,fng / Delta^2_N,sky avg modes:', 1. / signal_attenuation / n.sqrt(t_eff)

if False:
    p.subplot(122)
    p.plot(fng_bins * 1e3, fng_wgts)
    p.plot(n.array([fng_mn,fng_mn]) * 1e3, [0,1], 'k-')
    p.plot(n.array([fng_mx,fng_mx]) * 1e3, [0,1], 'k-')

p.show()
