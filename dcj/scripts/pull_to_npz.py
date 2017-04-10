
import aipy as a, numpy as n, pylab as p, sys, optparse
o = optparse.OptionParser()
o.set_usage('pull_to_npz.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,ant=True, pol=True, chan=True)
o.add_option('--accumulate',action='store_true',
    help="Concatenate all data into a single file instead of making individual files")

opts, args = o.parse_args(sys.argv[1:])
chan = int(opts.chan)
data = []
for uvfile in args:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    # Only select data that is needed to plot
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Read data from a single UV file
    for (uvw,t,(i,j)),d in uv.all():
        data.append(d[chan])
    if not opts.accumulate:
        print "saving to",uvfile+'.npz'
        n.save(uvfile+'.npz',n.array(data))
if opts.accumulate:
    n.save('data_slice.npz',n.array(data))
