#! /usr/bin/env python
import aipy as a, numpy as n, ephem
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True, src=True)
o.add_option('-P','--perceived', dest='perceived', action='store_true',
    help='Apply threshold for perceived Jys (using primary beam info)')
o.add_option('-S', '--standoff', dest='standoff', type='float', default=1.,
    help='Angle (in degrees) around a source to exclude locating lesser sources')
o.add_option('-t', '--thresh', dest='thresh', type='float', default=20.,
    help='Lower threshold in Jys for selecting a source.')
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa.set_jultime(2455746.)
srclist,cutoff,catalog = a.scripting.parse_srcs(opts.src, opts.cat)
cats = {}
for cat in catalog:
    cats[cat] = a.cal.get_catalog(opts.cal, srclist, cutoff, [cat])
    cats[cat].compute(aa)
decs = n.arange(-n.pi/2, n.pi/2, .001)
z = n.cos(decs - aa.lat)
y = n.sin(decs - aa.lat)
x = n.zeros_like(y)
if opts.perceived:
    print 'Computing primary beam response'
    bm = aa[0].bm_response((x,y,z), pol=opts.pol[0])
    bm *= aa[0].bm_response((x,y,z), pol=opts.pol[1])
    bm = bm.squeeze()
    bm_poly = n.polyfit(decs, n.log10(bm), deg=12)
    valid = n.where(n.abs(decs - aa.lat) <= n.pi/2, 1., 0)
#import pylab
#pylab.semilogy(decs, bm * valid, ',')
#pylab.semilogy(decs, 10**n.polyval(bm_poly, decs) * valid)
#pylab.show()

m = a.map.Map(fromfits=args[-1])
m.reset_wgt()
ras,decs = a.coord.eq2radec(m.px2crd(n.arange(m.map.map.size), ncrd=3))
if opts.perceived:
    print 'Applying beam response to map'
    valid = n.where(n.abs(decs - aa.lat) <= n.pi/2, 1., 0)
    gain = 10**n.polyval(bm_poly, decs) * valid
    m.map.map *= gain

print 'Selecting pixels'
srcmap = n.zeros_like(m.map.map)
px = n.where(m.map.map > opts.thresh)
jys = m.map.map[px]
ind = n.argsort(jys)[::-1]
jys = jys[ind]
px = px[0][ind]
print jys
print len(px)
ras,decs = a.coord.eq2radec(m.px2crd(px, ncrd=3))
prev_srcs = {}
colors = {}
for p,ra,dec in zip(px,ras,decs):
    src = a.amp.RadioFixedBody(ra=ra,dec=dec,jys=m.map.map[p],name=str(p))
    src.compute(aa)
    # Only select sources w/in 20-deg of zenith
    if n.abs(src.dec - aa.lat) > 20 * ephem.degree: continue
    flag = False
    for ep,esrc in prev_srcs.iteritems():
        esrc = esrc[0]
        if ephem.separation(src, esrc) < ephem.degree * opts.standoff:
            #print 'Source in px=%d is too close to %s' % (p, esrc.src_name)
            #print src._ra, src._dec, esrc._ra, esrc._dec
            prev_srcs[ep].append(src)
            srcmap[p] = colors[ep]
            flag = True
            break
    if flag: continue
    prev_srcs[p] = [src]
    colors[p] = n.random.randint(1,10)
    srcmap[p] = colors[p]
for p in px:
    try:
        srcs = prev_srcs[p]
        #if len(srcs) > 2: continue
        #print p, srcs[0]._ra, srcs[0]._dec, [src._jys for src in srcs]
        jys = [src._jys for src in srcs]
        # Exclude sources within 10-deg of galactic plane
        g = a.ephem.Galactic(srcs[0])
        if n.abs(g.lat) < 10 * ephem.degree and n.abs(g.lat) < 90 * ephem.degree: continue
        src_names = []
        src_sep = []
        inhelm = False
        for cat in cats:
          for src in cats[cat].values():
            try:
                #if ephem.separation(src,srcs[0]) < 20*ephem.arcminute:
                if ephem.separation(src,srcs[0]) < ephem.degree * opts.standoff:
                    if cat == 'helm': inhelm = True
                    #src_names.append(cat+'/'+src.src_name)
                    src_names.append(src.src_name)
                    src_sep.append(ephem.separation(src,srcs[0]))
            except(TypeError): pass
        index = n.argsort(src_sep)
        src_names = n.array(src_names)[index]
        src_sep = n.array(src_sep)[index]
        if False:
            print srcs[0]._ra, srcs[0]._dec,'|',
            print g.lat, max(jys), sum(jys)
            if inhelm: print '*',
            print '   ' + ' '.join([s+'/%04.1f'%(sep/ephem.arcminute) for s,sep in zip(src_names,src_sep)])
        else:
            if len(src_names) > 0:
                print src_names[0]
    except(KeyError): pass

m.map.map = srcmap
print 'Writing source map to srcmap.fits'
m.to_fits('srcmap.fits', clobber=True)
#print ','.join(['%s_%s' % (prev_srcs[p][0]._ra, prev_srcs[p][0]._dec) for p in px])


