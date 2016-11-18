#!/usr/bin/env python
print 'Hello World!'
import sys, capo as C, matplotlib.pyplot as plt, optparse, aipy as a
o = optparse.OptionParser()
o.set_usage('do_dataview.py [options] *.uv')
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('--max',default=0, dest='max',help='maximum value of colorbar')
o.add_option('--drng',default=4,dest='drng',help='dynamic range of colorbar')
o.add_option('--clean',default=1e-3,dest='clean',help='fractional tolerance for d/dr transforms')
o.add_option('--savepath',default=None,dest='save',help='If you want to save the figure, where? Give full path and name.png')
o.add_option('-v','--verbose',action='store_true',dest='verb',help='Toggle verbosity')
opts, args = o.parse_args(sys.argv[1:])
print args
mx,drng,clean = float(opts.max),float(opts.drng),float(opts.clean)

def bl2tup(bl):
    a1,a2 = map(int,bl.split('_'))
    return (a1,a2)
print 'Loading '
t,d,f = C.arp.get_dict_of_uv_data(args,antstr=opts.ant,polstr=opts.pol)

dd,ff = d[bl2tup(opts.ant)][opts.pol],f[bl2tup(opts.ant)][opts.pol]
print 'Plotting'
if not opts.save is None: C.sak.full_data_view(dd,ff,mx=mx,drng=drng,save=opts.save,clean=clean,verbose=opts.verb)
else: C.sak.full_data_view(dd,ff,mx=mx,drng=drng,clean=clean,verbose=opts.verb)

plt.show()
