#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import optparse, sys, math

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True, src=True)
o.add_option('-a', '--ant', dest='ant', type=int, default=0,
    help='Which antenna to use in the plot. Default = 0.')
o.add_option('-n', '--nbins', dest='nbins', type=float, default=96.,
    help='Number of points to evaluate beam response (over 24 hours).')
o.add_option('-b', '--binsize', dest='binsize', default='15',
    help='Number of minutes to include in each LST bin.  Default is 15.')
opts, args = o.parse_args(sys.argv[1:])

srclist = (opts.src).split(',')

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], .150, 1)
cat = a.cal.get_catalog(opts.cal, srclist)
cat.compute(aa)

twopi = 2 * math.pi
step = twopi/((60/float(opts.binsize))*24.)
bins = n.arange(0,twopi,step)
print len(bins)
timelist = (bins/twopi)+2455039
#timelist = n.arange(0,1,1/opts.nbins)+2455039
#timelist = n.arange(.91,.92,.0001)+2455039

bm = {}
lstlist = []
for k in cat:
    uv.rewind()
    print k
    bm[k] = []
    for time in timelist:
            aa.set_jultime(time)
            #lstlist.append(float(aa.sidereal_time()))
            cat.compute(aa)
            bm[k].append(aa[opts.ant].bm_response(cat[k].get_crds('top'), pol=opts.pol)[0][0])

for k in cat:
    print 'plotting', k
    p.plot(bins,cat[k].jys*(n.array(bm[k])**2),'.')
    #p.plot(timelist,n.sqrt(cat[k].jys)*bm[k],'.')

p.show()
