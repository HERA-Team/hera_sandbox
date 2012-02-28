#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a, sys, optparse, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o, max=True, drng=True)
opts,args = o.parse_args(sys.argv[1:])

plotdim1 = int(math.sqrt(len(args)))
plotdim2 = int(math.ceil(len(args)/plotdim1))
for cnt, filename in enumerate(args):
    srcname = filename.split('.')[-2].split('_')[-1]
    dfile = n.load(filename)
    d = n.log10(n.abs(dfile['srcest']))
    #d = n.log10(n.real(dfile['srcest']))
    #d = n.angle(dfile['srcest'])
    #d = n.abs(dfile['srcest'])
    if not opts.max is None: dmax = opts.max
    else: dmax = d.max()
    if not opts.drng is None: dmin = dmax - opts.drng
    else: dmin = d.min()
    p.subplot(plotdim1, plotdim2, cnt+1)
    p.title(srcname)
    p.imshow(d, aspect='auto', vmax=dmax, vmin=dmin)
    p.colorbar(shrink=.6)
p.show()
