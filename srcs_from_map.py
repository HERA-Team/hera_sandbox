#! /usr/bin/env python
import aipy as a, numpy as n, ephem
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-S', '--standoff', dest='standoff', type='float', default=1.,
    help='Angle (in degrees) around a source to exclude locating lesser sources')
o.add_option('-t', '--thresh', dest='thresh', type='float', default=20.,
    help='Lower threshold in perceived Jys for selecting a source.')
opts, args = o.parse_args(sys.argv[1:])

print 'Computing primary beam response'
aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa.set_jultime(2455016.)
decs = n.arange(-n.pi/2, n.pi/2, .001)
z = n.cos(decs - aa.lat)
y = n.sin(decs - aa.lat)
x = n.zeros_like(y)
bm = aa[0].bm_response((x,y,z), pol=opts.pol[0])
bm *= aa[0].bm_response((x,y,z), pol=opts.pol[1])
bm = bm.squeeze()
bm_poly = n.polyfit(decs, n.log10(bm), deg=12)
valid = n.where(n.abs(decs - aa.lat) <= n.pi/2, 1., 0)
#import pylab
#pylab.semilogy(decs, bm * valid, ',')
#pylab.semilogy(decs, 10**n.polyval(bm_poly, decs) * valid)
#pylab.show()

print 'Applying beam response to map'
m = a.map.Map(fromfits=args[-1])
m.reset_wgt()
ras,decs = a.coord.eq2radec(m.px2crd(n.arange(m.map.map.size), ncrd=3))
valid = n.where(n.abs(decs - aa.lat) <= n.pi/2, 1., 0)
gain = 10**n.polyval(bm_poly, decs) * valid
m.map.map *= gain

print 'Selecting pixels'
px = n.where(m.map.map > opts.thresh)
jys = m.map.map[px]
order = list(n.argsort(jys))
order.reverse()
order = n.array(order)
px = px[0].take(order)
jys = jys.take(order)
print px
print len(px)
ras,decs = a.coord.eq2radec(m.px2crd(px, ncrd=3))
prev_srcs = {}
for p,ra,dec in zip(px,ras,decs):
    src = a.amp.RadioFixedBody(ra=ra,dec=dec,jys=m.map.map[p],name=str(p))
    src.compute(aa)
    flag = False
    for ep,esrc in prev_srcs.iteritems():
        esrc = esrc[0]
        if ephem.separation(src, esrc) < ephem.degree * opts.standoff:
            #print 'Source in px=%d is too close to %s' % (p, esrc.src_name)
            #print src._ra, src._dec, esrc._ra, esrc._dec
            prev_srcs[ep].append(src)
            flag = True
            break
    if flag: continue
    prev_srcs[p] = [src]
for p in px:
    try:
        srcs = prev_srcs[p]
        print p, srcs[0]._ra, srcs[0]._dec, [src._jys for src in srcs]
    except(KeyError): pass

print ','.join(['%s_%s' % (prev_srcs[p][0]._ra, prev_srcs[p][0]._dec) for p in prev_srcs])


