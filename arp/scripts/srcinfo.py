#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys, ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--sep',dest='sep', type='float', default=1.,
    help='Include areas within the specified angular separation (in degrees) of any sources listed in --src.')
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.src.get_catalog(srclist, cutoff, catalogs)

for srcname in cat:
    src = cat[srcname]
    ephem.FixedBody.compute(src, ephem.J2000)
    
cats = {}
for catname in catalogs:
    cats[catname] = a.src.get_catalog(catalogs=[catname])
    for srcname in cats[catname]:
        try: ephem.FixedBody.compute(cats[catname][srcname], ephem.J2000)
        except(TypeError): pass

srcmeas = {}
for filename in args:
    data = [L.split() for L in open(filename).readlines()]
    data = [L[:1] + map(float, L[1:]) for L in data if len(L) == 5]
    for d in data:
        try: srcmeas[d[0]].append((d[1], d[4]))
        except(KeyError): srcmeas[d[0]] = [(d[1], d[4])]
    

fqs = 10**n.arange(-2, 0, .01)
for s1name in cat:
    s1 = cat[s1name]
    s1x,s1y,names = [], [], []
    for cat2name in cats:
        cat2 = cats[cat2name]
        for s2name in cat2:
            s2 = cat2[s2name]
            try:
                if ephem.separation(s1, s2) <= opts.sep * a.img.deg2rad:
                    s1x.append(s2.mfreq)
                    s1y.append(s2._jys)
                    if s2.index != 0:
                        s2.update_jys(fqs)
                        p.loglog(fqs, s2.jys, ':', label=(s2name+'(%s)'%cat2name))
                    names.append(s2name + '(%s)' % cat2name)
                    if srcmeas.has_key(s2name):
                        data = n.array(srcmeas[s2name])
                        p.loglog(data[:,0], data[:,1], 'k^', label=s2name+'(meas)')
            except(TypeError): pass
    if srcmeas.has_key(s1name):
        data = n.array(srcmeas[s1name])
        p.loglog(data[:,0], data[:,1], 'k^', label=s1name+'(meas)')
    p.loglog(s1x, s1y, 'k.')
    p.title(s1name)
    p.grid()
    p.legend()
    p.xlim(.050, .250)
    for x,y,n in zip(s1x,s1y,names): p.text(x, y, n)
    p.show()
  
    
