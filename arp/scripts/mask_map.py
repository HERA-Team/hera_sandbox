#! /usr/bin/env python
import aipy as a, numpy as n, ephem
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa.set_jultime(2455747.)
cen = a.cal.get_catalog(opts.cal, ['cen'], catalogs=['misc']).values()[0]
cen.compute(aa)
cra,cdec = cen.get_crds('eq', ncrd=2)
for ifile in args:
    print 'Reading', ifile
    m = a.map.Map(fromfits=ifile)

    print 'Generating pixel coordinates'
    exyz = m.px2crd(n.arange(m.map.map.size), ncrd=3)
    eq2ga = a.coord.convert_m('eq','ga')
    gxyz = n.dot(eq2ga, exyz)
    ras,decs = a.coord.eq2radec(exyz)
    glon,glat = a.coord.eq2radec(gxyz)

    print 'Applying mask'
    valid = n.where(n.abs(decs - aa.lat) < 20 * a.ephem.degree, 1, 0)
    valid *= n.where(n.abs(glat) < 10 * a.ephem.degree, 0, 1)
    #ex,ey,ez = exyz
    #cx,cy,cz = cen_xyz
    #r = n.sqrt((ex-cx)**2 + (ey-cy)**2 + (ez-cz)**2)
    #valid *= n.where(r < 5 * a.ephem.degree, 0, 1)
    #valid *= n.where(r < 5 * a.ephem.degree, 0, 1)
    valid *= n.logical_or(n.where(n.abs(ras-cra) < 2.5 * a.ephem.degree,0,1), n.where(n.abs(decs-cdec) < 5 * a.ephem.degree, 0, 1))
    m.map.map *= valid
    m.wgt.map *= valid

    ofile = 'mask_' + ifile
    print 'Writing masked map to', ofile
    m.to_fits(ofile, clobber=True)


