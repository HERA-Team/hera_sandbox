#! /usr/bin/env python
import aipy as a, numpy as n, pylab as P
import capo as C
import sys, optparse
#from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
NCHAN = uv['nchan']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

STEP = 37
us,ds = {}, {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    times = []
    uv.select('auto', -1, -1, include=False)
    uv.select('antennae', 40, -1, include=False)
    uv.select('antennae', 55, -1, include=False)
    #uv.select('decimate', 100, 10)
    uv.select('decimate', 16, 2)
    for (uvw, t, (i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            print t, len(times)
        #print i,j
        u,v,w = aa.gen_uvw(i,j, src='z')
        umag = n.sqrt(u**2 + v**2).squeeze()
        d /= aa.passband(i,j)
        for CH in range(0, NCHAN, STEP):
            if not us.has_key(CH): us[CH], ds[CH] = {}, {}
            valid = n.where(d[CH:CH+STEP] == 0, 0, 1)
            for u,_d in zip(umag[CH:CH+STEP].compress(valid), n.abs(d[CH:CH+STEP].compress(valid))):
                bin = n.around(n.log10(u), 1)
                if not us[CH].has_key(bin): us[CH][bin], ds[CH][bin] = [],[]
                us[CH][bin].append(u)
                ds[CH][bin].append(_d)

for color, CH in zip('kgbrcmy', range(0, NCHAN, STEP)):
    bins = ds[CH].keys()
    #for b in bins:
    #    us[CH][b] = n.concatenate(us[CH][b])
    #    ds[CH][b] = n.concatenate(ds[CH][b])
    stds = n.array([n.std(ds[CH][b]) / n.sqrt(len(ds[CH][b])) for b in bins])
    values = n.array([n.average(ds[CH][b]) for b in bins])
    P.errorbar(10**n.array(bins), values, 2*stds, fmt=color+'.')
ax = P.gca(); ax.set_xscale('log'); ax.set_yscale('log')
P.show()
