#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re

re_filename = re.compile(r'.*_c(\d+)_(\d+)_.*')

o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True)
o.add_option('--sfreq', dest='sfreq', type='float', default=0.1212890625)
o.add_option('--sdf', dest='sdf', type='float', default=0.00029296875)
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff = a.scripting.parse_srcs(opts.src)
cat = a.src.get_catalog(srclist, cutoff)

def ch2freq(chan, sfreq=0.1212890625, sdf=0.00029296875):
    return (sfreq + chan*sdf)

srccrd = {}
fqs = []
srcspec = {}
for src in cat:
    srccrd[src] = a.coord.radec2eq((cat[src]._ra, cat[src]._dec))
    srcspec[src] = []

for map in args:
    print 'Reading', map,
    ch1, ch2 = re_filename.match(map).groups()
    ch = .5 * (float(ch1) + float(ch2))
    fq = ch2freq(ch, sfreq=opts.sfreq, sdf=opts.sdf)
    fqs.append(fq)
    print 'at %0.4f GHz' % fq
    m = a.map.Map(fromfits=map)
    if not opts.nside is None:
        nm = a.healpix.HealpixMap(nside=opts.nside)
        nm.from_hpm(m)
        m = nm
    m.set_interpol(False)
    for src in cat:
        srcspec[src].append(m[tuple(srccrd[src])])

import pylab
fqs = n.array(fqs)
inds = n.argsort(fqs)
fqs = fqs[inds]
for src in cat:
    #print "'%s' : {'jys':10**%f, 'index':[%f], }," % (src
    srcspec[src] = n.array(srcspec[src])[inds]
    pylab.plot(fqs, srcspec[src], label=src)
pylab.legend()
pylab.show()
