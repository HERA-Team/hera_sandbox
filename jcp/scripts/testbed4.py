#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import ephem, sys, optparse, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('-c', '--chan', dest='chan', type='int', default=600,
    help='Channel')
o.add_option('-d', '--decimate', dest='decimate', default=1,
    help='Decimate the actual data')
opts,args = o.parse_args(sys.argv[1:])

aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9'))

auto = []
lsts = []
my_auto = None
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    uv.select('antennae', 1, 1)
    uv.select('decimate', opts.decimate, 0)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        auto.append(d[opts.chan])
        aa.set_jultime(t)
        print aa.sidereal_time() - uv['lst']
        lsts.append(aa.sidereal_time())
        #sys.exit()

lsts=n.array(lsts)

p.plot((24*lsts)/(2*math.pi), auto)
p.show()
