#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, math

o = optparse.OptionParser()
o.set_usage('plot_topo_hmap.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, chan=True, cmap=True, max=True, drng=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
opts, args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap) 
m2 = int(math.sqrt(len(args)))
m1 = int(math.ceil(float(len(args)) / m2))

for cnt, filename in enumerate(args):
    print 'Reading',filename
    h = a.map.Map(fromfits=filename)
    p.subplot(m2, m1, cnt+1)
    C.jcp.plot_hpx(h, mode=opts.mode, mx=opts.max, drng=opts.drng, colorbar=True)
    p.title(filename)
p.show()
