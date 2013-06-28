#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p, optparse, sys, os

o = optparse.OptionParser()
o.set_usage('plot_srctrack.py npzfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True)

opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))

srclist = []
for file in args:
    srclist.append(os.path.basename(file).split('__')[0])
srclist = set(srclist)

nplots = len(srclist)
ncols = n.ceil(n.sqrt(nplots))
nrows = n.round(n.sqrt(nplots))


for i, src in enumerate(srclist):
    #p.subplot(nrows,ncols,i+1)
    p.title(src)
    files = [s for s in args if src in s]
    for file in files:
        print file
        npz = n.load(file)
        times = []
        for t in npz['times']:
            aa.set_jultime(t)
            times.append(aa.sidereal_time())
        spec = n.abs(npz['spec'])
        p.semilogy(times,spec,'-', label=src)
p.legend()
p.show()
    
