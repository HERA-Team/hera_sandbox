#! /usr/bin/env python
import numpy as n
import aipy as a
import sys, os, optparse 

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--scale', type=float, 
            help='multiply the data by this number')
opts, args = o.parse_args(sys.argv[1:])
if opts.scale is None:
   print "ERROR: no scale supplied."
   sys.exit()


def mfunc(uv, p, d, f):
    uvw,t,(ij) = p
    d = d*opts.scale
    return p, d, f
     

for infile in args:
    outfile = infile+'S'
    print infile, '-->', outfile
    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue
    
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc, raw=True, append2hist='SCALE_DATA: ' + ' '.join(sys.argv) + '\n')
