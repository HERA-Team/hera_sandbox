#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

exec('import %s as cal' % opts.cal)
afreqs = n.arange(.130, .190, .01)
srclist = cal.src_prms.keys() + ['cyg']

cat = a.src.get_catalog(srclist)
try: del(cat['Sun'])
except(KeyError): pass
#del(cat['hyd'])
srclist = cat.keys()
data = {}
good_list = []
for src in srclist:
    cat[src].update_jys(afreqs)
    flx1 = cat[src].jys
    if n.all(flx1 < 1e-2): continue
    good_list.append(src)
    try: cat[src].set_params(cal.src_prms[src])
    except(KeyError): pass
    cat[src].update_jys(afreqs)
    flx2 = cat[src].jys
    data[src] = (cat[src]._ra, cat[src]._dec, flx1, flx2)
    print src, n.average(flx2/flx1)

srclist = good_list

ras = n.array([data[s][0] for s in srclist])
inds = n.argsort(ras)
ras = ras[inds]
decs = n.array([data[s][1] for s in srclist])[inds]
flx1s = n.array([data[s][2] for s in srclist])[inds]
flx2s = n.array([data[s][3] for s in srclist])[inds]

for i, af in enumerate(afreqs):
    p.subplot(2, 3, i+1)
    p.loglog(flx1s[:,i], flx2s[:,i], 'k.')
    x = 10**n.arange(0,5,.01)
    f = n.average(flx2s[:,i]/flx1s[:,i])
    p.loglog(x, f*x, 'k:', label='f=%1.3f' % (f))
    p.loglog(x, x, 'k-')
    p.title('%3.1f MHz' % (1e3*af))
    p.xlabel('Catalog Flux Density (Jy)')
    p.ylabel('Measured Flux Density (Jy)')
    p.xlim(1,1e5)
    p.ylim(1,1e5)
    p.grid()
    p.legend(loc='best')

p.show()
'''

p.subplot(121)
for i, af in enumerate(afreqs):
    p.plot(ras * 12 / n.pi , flx2s[:,i]/flx1s[:,i], '.')
p.xlabel('Right Ascension (hours)')
p.ylabel('Measured Flux / Catalog Flux')
p.xlim(0,24)

p.subplot(122)
for i, af in enumerate(afreqs):
    p.plot(decs * 180 / n.pi , flx2s[:,i]/flx1s[:,i], '.')
p.xlabel('Declination (degrees)')
p.ylabel('Measured Flux / Catalog Flux')
p.xlim(0,90)

p.show()
'''
