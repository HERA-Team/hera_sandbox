#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True, cmap=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, n.array([.1375]))

ants = {}
times = {}
for filename in args:
    print 'Reading', filename
    src = filename.split('.')[-1]
    src = a.cal.get_catalog(opts.cal, [src]).values()[0]
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, 'auto', opts.pol)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if not times.has_key(t):
            aa.set_jultime(t)
            src.compute(aa)
            if src.alt < .5: continue
            times[t] = {'top':src.get_crds('top', ncrd=3)}
        times[t][i] = d[0]
        ants[i] = None

tlist = times.keys()
for t in tlist:
    for i in ants:
        if not times[t].has_key(i):
            del(times[t])
            break
tlist = times.keys()
tlist.sort()
for i in ants: ants[i] = n.array([times[t][i] for t in tlist])
top = n.array([times[t]['top'] for t in tlist])
tx,ty,tz = top.transpose()

a0 = sum([n.abs(ants[i]) for i in ants]) / len(ants)
a0g = a0.sum()
for i in ants:
    ai = n.abs(ants[i])
    gain = ai.sum() / a0g
    ai /= a0 * gain
    p.semilogy(ai.clip(.25,4))
p.show()
