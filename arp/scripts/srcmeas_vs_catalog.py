#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys, ephem

o = optparse.OptionParser()
o.add_option('--cat', dest='cat', default='helm,misc',
    help='A list of catalogs to use.')
o.add_option('--xcat', dest='xcat', default='',
    help='A list of included catalogs to exclude from plotting.')
o.add_option('--sep',dest='sep', type='float', default=1.,
    help='Include areas within the specified angular separation (in degrees) of any sources listed in --src.')
o.add_option('--sum',dest='sum', action='store_true',
    help='Sum all sources within "sep" of the position of the measured source.  Otherwise, take brightest source.')
o.add_option('--ratio',dest='ratio', action='store_true',
    help='Plot ratio of measured flux to catag flux.  Otherwise, make a flux-flux plot.')
opts,args = o.parse_args(sys.argv[1:])

catalogs = opts.cat.split(',')
xcatalogs = opts.xcat.split(',')
pcats = [c for c in catalogs if not c in xcatalogs]

srcmeas = {}
for filename in args:
    data = [L.split() for L in open(filename).readlines()]
    data = [L[:1] + map(float, L[1:]) for L in data if len(L) == 5]
    for d in data:
        try: srcmeas[d[0]].append((d[1], d[4]))
        except(KeyError): srcmeas[d[0]] = [(d[1], d[4])]

srclist = srcmeas.keys()
cat = a.src.get_catalog(srclist, catalogs=catalogs)
for srcname in cat:
    src = cat[srcname]
    ephem.FixedBody.compute(src, ephem.J2000)

cats = {}
for catname in catalogs:
    cats[catname] = a.src.get_catalog(catalogs=[catname])
    for srcname in cats[catname]:
        try: ephem.FixedBody.compute(cats[catname][srcname], ephem.J2000)
        except(TypeError): pass

equal = 10**(n.arange(0,5,.01))
mfreqs = n.array([.14,.15,.16,.17])
sub2 = int(n.sqrt(len(pcats)))
sub1 = int(n.ceil(len(pcats) / float(sub2)))
badsrcs = {}
for i,cat2name in enumerate(pcats):
    p.subplot(sub1,sub2, i+1)
    cat2 = cats[cat2name]
    catflx, names = [], []
    mfreq = None
    for s1name in cat:
        s1 = cat[s1name]
        maxsrc = 0
        for s2name in cat2:
            s2 = cat2[s2name]
            try:
                if ephem.separation(s1, s2) <= opts.sep * a.img.deg2rad:
                    if s2._jys > maxsrc:
                        if opts.sum: maxsrc += s2._jys
                        else: maxsrc = s2._jys
            except(TypeError): pass
        if maxsrc > 0:
            catflx.append(maxsrc)
            names.append(s1name)
            mfreq = s2.mfreq
    j = n.argmin(n.abs(mfreqs - mfreq))
    catflx = n.array(catflx)
    measflx = n.array([srcmeas[s][j][1] for s in names])
    badsrc = n.argwhere(n.abs(measflx/catflx - 1) > .25)
    for b in badsrc:
        b = names[b]
        if not badsrcs.has_key(b): badsrcs[b] = []
        badsrcs[b].append(cat2name)
    if opts.ratio:
        p.semilogx(catflx, measflx/catflx, 'k.')
        p.ylabel('meas/cat (%1.3f)' % mfreqs[j])
        #for b in badsrc: p.text(catflx[b], measflx[b]/catflx[b], names[b])
        #p.ylim(.5,1.5)
    else:
        p.loglog(catflx, measflx, 'k.')
        p.loglog(equal,equal, 'k:')
        p.loglog(equal,1.25*equal, 'r:')
        p.loglog(1.25*equal,equal, 'r:')
        p.ylabel('meas (%1.3f)' % mfreqs[j])
        #for b in badsrc: p.text(catflx[b], measflx[b], names[b])
        p.ylim(2,2e4)
    p.xlabel(cat2name+' (%1.3f)' % mfreq)
    p.grid()
    p.xlim(2,2e4)
for b in badsrcs:
    if len(badsrcs[b]) > 1: print b, badsrcs[b]
p.show()
