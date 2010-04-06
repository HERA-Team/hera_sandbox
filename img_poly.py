#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, loc=True)
o.add_option('--sfreq', dest='sfreq', type='float', default=0.1212890625,
    help='Start frequency of channels (GHz).')
o.add_option('--sdf', dest='sdf', type='float', default=0.00029296875,
    help='Frequency interval between channels (GHz).')
o.add_option('-d', '--deg', dest='deg', type='float', default=4,
    help='Degree of polynomial to use.')
opts,args = o.parse_args(sys.argv[1:])

re_ch = re.compile('.*c(\d+)_(\d+).*', re.IGNORECASE)
freqs = []
imgs = []
for filename in args:
    print 'Reading %s' % filename
    c1,c2 = map(int, re_ch.match(filename).groups())
    c = .5 * (c1 + c2)
    freqs.append(opts.sfreq + c * opts.sdf)
    print '   Freq:', freqs[-1]
    img,kwds = a.img.from_fits(filename)
    imgs.append(img.flatten().astype(n.float32))

freqs = n.array(freqs, dtype=n.float32)
imgs = n.array(imgs, dtype=n.float32)
print imgs.shape

print 'Fitting polynomials...'
print freqs.dtype, imgs.dtype
polys = n.polyfit(freqs, imgs, opts.deg)
print 'Done.'
simgs = []

print 'Smoothing imgs...'
for i,f in enumerate(freqs):
    print i
    f = f**n.arange(opts.deg, -1, -1)
    f.shape = (1,f.size)
    print polys.shape, f.shape
    x = n.dot(f, polys).squeeze()
    print x.shape
    simgs.append(x)
    imgs[i] -= x
simgs = n.array(simgs, n.float32)
print 'Done.'

print 'Saving...'
for i,filename in enumerate(args):
    print '    ', filename
    img,kwds = a.img.from_fits(filename)
    im = simgs[i] ; im.shape = img.shape
    a.img.to_fits(filename[:-8]+'sm'+filename[-8:], im, **kwds)
    im = imgs[i] ; im.shape = img.shape
    a.img.to_fits(filename[:-8]+'rs'+filename[-8:], im, **kwds)
print 'Done.'
