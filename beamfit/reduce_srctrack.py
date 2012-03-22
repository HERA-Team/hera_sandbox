#! /usr/bin/env python
import aipy as a, numpy as n, optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--cat', dest='cat', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
             
opts,args = o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal, n.array([.150]))

for filename in args:
    print 'Reading', filename
    src = os.path.basename(filename).split('__')[0]
    srclist, cutoff, catalogs = a.scripting.parse_srcs(src, opts.cat)
    src = a.cal.get_catalog(opts.cal,srclist)[src]
    npz = n.load(filename)
    x,y,z = [],[],[]
    for t in npz['times']:
        aa.set_jultime(t)
        src.compute(aa)
        xi,yi,zi = src.get_crds('top')
        x.append(xi)
        y.append(yi)
        z.append(zi)
    x,y,z = n.array(x), n.array(y), n.array(z)
    dict = {}
    for file in npz.files: dict[file] = npz[file]
    n.savez('_'+filename, x=x, y=y, z=z, **dict)

