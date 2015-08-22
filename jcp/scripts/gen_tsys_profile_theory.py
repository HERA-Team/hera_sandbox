#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import ephem, sys, optparse, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-m', '--map', dest='map',
    help='Haslam map; or other!.')
o.add_option('-c', '--chan', dest='chan', type='int', default=600,
    help='Channel')
o.add_option('--mfreq', dest='mfreq', type='float',
    help='Frequency of map in GHz.')
o.add_option('--spind', dest='spind', default=-2.52, type='float',
    help='Spectral index to use when extrapolating map to other frequencies.')
o.add_option('-P', '--power', dest='pflag', action='store_true',
    help='Put profile in units of total power.')
opts,args = o.parse_args(sys.argv[1:])


fq = 0.03
aa = a.cal.get_aa(opts.cal,n.array([fq]))
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

lsts = []
dat = []
for ha in n.arange(0,2*n.pi, .1):
    lsts.append(ha)
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
#time change correction?
lsts_fit = n.array(lsts + lsts + lsts)
lsts = n.array(lsts)
dat_fit = n.concatenate([dat, dat, dat])
sync_poly = n.polyfit(lsts_fit, dat_fit, deg=12)
sync_auto = n.polyval(sync_poly, lsts)

#output the data to .npz file
#n.savez("tsys_profile.npz",lsts=lsts, sync_auto=sync_auto)

taxis = (lsts * 24)/(n.pi*2)
p.plot(taxis, sync_auto,label = 'Model Sky')#tweaked to match plot_uv.py old:(*12 /n.pi)
p.xlabel("Time (Hours)")
p.ylabel("Temperature (K)")
p.show()
