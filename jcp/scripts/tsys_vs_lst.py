#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import ephem, sys, optparse, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-m', '--map', dest='map',
    help='Haslam map; or other!.')
o.add_option('--mfreq', dest='mfreq', type='float',
    help='Frequency of map in GHz.')
o.add_option('--spind', dest='spind', default=-2.52, type='float',
    help='Spectral index to use when extrapolating map to other frequencies.')
opts,args = o.parse_args(sys.argv[1:])

sdf = .1
fq = opts.mfreq
aa = a.cal.get_aa(opts.cal, sdf, fq, 1)
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
    print ha, aa.lat
    ex, ey, ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten() 
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = h[ex,ey,ez] * (fq/opts.mfreq)**opts.spind
    tsky.shape = resp.shape
    tsky = n.where(invalid, 0, tsky)
    dat.append(n.sum(tsky * resp) / n.sum(resp))

dat = n.array(dat)
lsts_fit = n.array(lsts_fit + lsts_fit + lsts_fit)
dat_fit = n.concatenate([dat, dat, dat])
#sync_poly = n.polyfit(lsts_fit, dat_fit, deg=12)
#sync_auto = n.polyval(sync_poly, lsts)

if False:
    poly = n.polyfit(sync_auto, auto, deg=1)
    print poly
    gain, T_rx = poly
    T_rx /= gain
    print gain, T_rx
#else:
    #gain, T_rx = 1057., 92.7
    #gain, T_rx = 1143., 76.7

#output the data to .npz file
#n.savez("tsys_profile.npz",lsts=lsts, sync_auto=sync_auto)

p.plot(n.arange(0,2*n.pi, .1) * 12 / n.pi, dat)
p.xlabel('LST')
p.ylabel('Tsky (K)')
p.show()
'''

#nplots=2
#if opts.pflag: nplots+=1
if opts.temp: nplots+=1
p.subplot(nplots,1,1)
taxis = (lsts * 24)/(n.pi*2)
p.plot(taxis, sync_auto,label = 'Model Sky')#tweaked to match plot_uv.py old:(*12 /n.pi)
p.plot(taxis, auto/gain - T_rx, '.',label='Total Power (Autocorrelation)')
p.legend(loc=2)
p.xlabel("Time (Hours)")
p.ylabel("Temperature (K)")
p.subplot(nplots,1,2)
#p.plot(lsts, auto, '.')
p.plot(sync_auto, auto/gain, '.')
Ts = n.arange(0,500,1)
p.plot(Ts, (Ts+T_rx), ':')
p.xlabel = ("Model Temperature")
p.ylabel = ("Autocorrelation Temperature")
if opts.temp: 
    p.subplot(nplots,1,3)
    p.plot(temp['lsts'],temp['T'])
p.show()
if opts.pflag:
    p.close('all')
    p.plot(lsts, power_auto)
    p.show()



resp = a.img.recenter(resp, (100,100))
tsky = a.img.recenter(tsky, (100,100))

print n.sum(tsky * resp) / n.sum(resp)

p.subplot(131)
p.imshow(resp)
p.colorbar(shrink=.5)

p.subplot(132)
p.imshow(tsky)
p.colorbar(shrink=.5)

p.subplot(133)
p.imshow(tsky * resp)
p.colorbar(shrink=.5)

p.show()
'''
