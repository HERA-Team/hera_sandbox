#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

uv = a.miriad.UV(sys.argv[1])
aa = a.cal.get_aa('pgb966', uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

dat = {}
wgt = {}

for filename in sys.argv[1:]:
    print filename
    srcname = filename.split('.')[-1][3:]
    cat = a.cal.get_catalog('pgb966', [srcname])
    src = cat[srcname]
    start_t = float('.'.join(filename.split('.')[1:3]))
    aa.set_jultime(start_t+.02) # halfway into file
    src.compute(aa)
    if src.alt < a.img.deg2rad * 30: continue
    uv = a.miriad.UV(filename)
    for (uvw,t,bl),d,f in uv.all(raw=True):
        aa.set_jultime(t)
        src.compute(aa)
        s_eq = cat.get_crds('eq', ncrd=3)
        aa.sim_cache(s_eq)
        w = aa.bm_response(0,0,pol='yy').squeeze()
        w *= n.sqrt(src.get_jys())
        w *= n.logical_not(f)
        w *= w
        rnd_t = n.round(t/4, 2) * 4
        dat[rnd_t] = dat.get(rnd_t,0) + (d * w).sum()
        wgt[rnd_t] = wgt.get(rnd_t,0) + w.sum()
    del(uv)
        
ts = dat.keys(); ts.sort()
ts = n.array(ts)
dats = n.array([dat[t] for t in ts])
wgts = n.array([wgt[t] for t in ts])
valid = n.where(wgts == 0, 0, 1)
dats = dats.compress(valid)
wgts = wgts.compress(valid)
ts = ts.compress(valid)

dats /= wgts
print repr(ts)
print repr(dats)
p.plot(ts, dats)
#p.ylim(0,2)
p.show()
