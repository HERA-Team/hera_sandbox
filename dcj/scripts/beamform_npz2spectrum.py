#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, scipy.optimize
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, pol=True)
o.add_option('-s', '--src', dest='src', type='str', 
    help='Source to use for calibration.')
o.add_option('--cat', dest='cat', type='str', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
opts,args = o.parse_args(sys.argv[1:])

def filename2src(f):
    #return f.split('_')[-1]
    return f.split('_')[-1].split('.')[0]

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.select_chans(chans)

if opts.src != None:
    srclist,cutoff,catalogs, = a.scripting.parse_srcs(opts.src, opts.cat)
    calsrc = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
    calsrc = calsrc.values()[0]
else:
    calsrc = None
srclist = [filename2src(f) for f in args]
print srclist
if not calsrc is None: assert(calsrc.src_name in srclist)
srclist,cutoff,catalogs, = a.scripting.parse_srcs(','.join(srclist), opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)

for srctrack in myfiles:
    print srctrack
    srcname = filename2src(srctrack)
    tracks[srcname] = n.load(srctrack)
