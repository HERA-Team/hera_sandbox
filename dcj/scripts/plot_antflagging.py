#! /usr/bin/env python
import aipy as a
import numpy as n
from pylab import *
import optparse, sys, os
from astropy.time import Time

o=optparse.OptionParser()
o.set_usage("plot_antflagging.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-v',action='store_true',help='turn on more verbs')
o.add_option('--plot_date',action='store_true',help='interpret column 1 as a jd and plot accordingly (otherwise \
ignore)')
opts,args=o.parse_args(sys.argv[1:])

ax = subplot(111)
for arg in args:
    D = n.loadtxt(arg)
    t = Time(D[:,0],scale='utc',format='jd')
    for i in n.arange(D.shape[1]-1):
        if opts.plot_date:
            plot_date(t.plot_date,D[:,i+1],'-',label=str(i))
        else:
            plot(D[:,i+1],'-',label=str(i))
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.8])
legend(loc='upper center',ncol=16,bbox_to_anchor=(0.5,-0.05),numpoints=1,fontsize='x-small')
figure()
show()
