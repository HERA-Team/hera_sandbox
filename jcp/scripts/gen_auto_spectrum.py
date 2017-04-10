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

ha = (20./24.)*2*n.pi
freqs = n.linspace(.1,.2,203)
dat = []
for fq in freqs:
    aa = a.cal.get_aa(opts.cal,fq)
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
    
    ex, ey, ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten() 
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = h[ex,ey,ez] * (fq/opts.mfreq)**opts.spind
    tsky.shape = resp.shape
    tsky = n.where(invalid, 0, tsky)
    dat.append(n.sum(tsky * resp) / n.sum(resp))
 
dat = n.array(dat)

if True:
    dat -= dat[102] * (freqs/freqs[102])**opts.spind
    #_dat = dat[102] * (freqs/freqs[102])**opts.spind
    #p.plot(freqs, _dat,label = 'Model Sky')

p.plot(freqs, dat,label = 'Model Sky')
p.xlabel("Frequency (GHz)")
p.ylabel("Temperature (K)")
p.show()
