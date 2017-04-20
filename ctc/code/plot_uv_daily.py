#! /usr/bin/env python

import numpy
import aipy, capo
import pylab as p
import glob
import optparse
import sys, os

###
# Plots UV files on a day-by-day basis
###

o = optparse.OptionParser()
o.set_usage('plot_uv_daily.py *npz')
aipy.scripting.add_standard_options(o,ant=True,pol=True)
o.add_option('--drng', default=3.0, type='float')
o.add_option('--mx', default=0.0, type='float')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

plotnum = 1
jd0 = args[0].split('/')[-1].split('.')[1]
dailyfiles = []
for file in args:
    jd = file.split('/')[-1].split('.')[1]
    if jd == jd0:
        dailyfiles.append(file)
    else:
        p.figure(plotnum)
        t,d,f = capo.miriad.read_files(dailyfiles,antstr=opts.ant,polstr=opts.pol,verbose=True)
        bl = (int(opts.ant.split('_')[0]),int(opts.ant.split('_')[1]))
        capo.plot.waterfall(d[bl][opts.pol],drng=opts.drng,mx=opts.mx,extent=(0,202,t['lsts'][-1]*12/numpy.pi,t['lsts'][0]*12/numpy.pi))
        p.title(jd0)
        p.xlabel('Channel')
        p.ylabel('LST Hour')
        jd0 = jd
        dailyfiles = []
        plotnum += 1
        dailyfiles.append(file)
p.figure(plotnum)
t,d,f = capo.miriad.read_files(dailyfiles,antstr=opts.ant,polstr=opts.pol,verbose=True)
bl = (int(opts.ant.split('_')[0]),int(opts.ant.split('_')[1]))
capo.plot.waterfall(d[bl][opts.pol],drng=opts.drng,mx=opts.mx,extent=(0,202,t['lsts'][-1]*12/numpy.pi,t['lsts'][0]*12/numpy.pi))
p.title(jd0)
p.xlabel('Channel')
p.ylabel('LST Hour')
p.show()
