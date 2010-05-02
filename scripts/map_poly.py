#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re

o = optparse.OptionParser()
o.add_option('--sfreq', dest='sfreq', type='float', default=0.1212890625,
    help='Start frequency of channels (GHz).')
o.add_option('--sdf', dest='sdf', type='float', default=0.00029296875,
    help='Frequency interval between channels (GHz).')
o.add_option('-d', '--deg', dest='deg', type='float', default=4,
    help='Degree of polynomial to use.')
opts,args = o.parse_args(sys.argv[1:])

re_ch = re.compile('.*c(\d+)_(\d+).*', re.IGNORECASE)
freqs = []
maps = []
for filename in args:
    print 'Reading %s' % filename
    c1,c2 = map(int, re_ch.match(filename).groups())
    c = .5 * (c1 + c2)
    freqs.append(opts.sfreq + c * opts.sdf)
    print '   Freq:', freqs[-1]
    m = a.map.Map(fromfits=filename)
    m.reset_wgt()
    maps.append(m.map.map.astype(n.float32))

freqs = n.array(freqs, dtype=n.float32)
maps = n.array(maps, dtype=n.float32)
print maps.shape

print 'Fitting polynomials...'
print freqs.dtype, maps.dtype
polys = n.polyfit(freqs, maps, opts.deg)
print 'Done.'
smaps = []

print 'Smoothing maps...'
for i,f in enumerate(freqs):
    print i
    f = f**n.arange(opts.deg, -1, -1)
    f.shape = (1,f.size)
    print polys.shape, f.shape
    x = n.dot(f, polys).squeeze()
    print x.shape
    smaps.append(x)
    maps[i] -= x
smaps = n.array(smaps, n.float32)
print 'Done.'

print 'Saving...'
for i,filename in enumerate(args):
    print '    ', filename
    m.map.map = smaps[i]; m.to_fits('sm_'+filename)
    m.map.map = maps[i]; m.to_fits('rs_'+filename)
print 'Done.'
