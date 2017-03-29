
import aipy as a, numpy as n, pylab as p, sys, optparse
import capo
o = optparse.OptionParser()
o.set_usage('uv_meanabs.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,ant=True, pol=True, chan=True,dec=True)
#o.add_option('--accumulate',action='store_true',
#    help="Concatenate all data into a single file instead of making individual files")

opts, args = o.parse_args(sys.argv[1:])
info, data, flag = capo.miriad.read_files(args,opts.ant,opts.pol,decimate=opts.decimate,recast_as_array=True)
D = []
for bl in data:
    for pol in data[bl]:
        D.append(data[bl][pol])
print n.mean(n.abs(D))

