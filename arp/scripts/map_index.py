#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re

o = optparse.OptionParser()
o.add_option('--sfreq', dest='sfreq', type='float', default=0.1001953125,
    help='Start frequency of channels (GHz).')
o.add_option('--sdf', dest='sdf', type='float', default=0.000390625,
    help='Frequency interval between channels (GHz).')
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
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
    if not opts.nside is None:
        nm = a.map.Map(nside=opts.nside)
        nm.from_map(m)
        m = nm
    m.reset_wgt()
    maps.append(m.map.map.astype(n.float32))

freqs = n.array(freqs, dtype=n.float32)
maps = n.array(maps, dtype=n.float32)
print maps.shape
avg = n.average(maps, axis=0)

print 'Fitting index...'
print freqs.dtype, maps.dtype
poly = n.polyfit(n.log10(freqs/.150), n.log10(n.abs(maps).clip(1e-6,n.Inf)), 1)
#poly = n.polyfit(freqs - .150, maps, 1)
print poly.shape
print 'Done.'

print 'Saving...'
#m.map.map = n.where(n.abs(poly[1]) > 1, poly[0]/poly[1], 0)
m.map.map = n.where(n.abs(avg) > 5, poly[0], 0)
m.to_fits('index_'+filename, clobber=True)
