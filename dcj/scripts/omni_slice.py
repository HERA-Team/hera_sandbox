#! /usr/bin/env python

import omnical
import aipy as a
import pylab
import numpy as n
import capo
import pickle
import matplotlib.pyplot as plt
from pylab import *
import optparse
import os, sys
from capo.dcj import file2jd
from capo.omni import from_npz,to_npz
### Functions ###
def dB(x):
    return  10*n.ma.log10(n.ma.abs(x))
### Options ###
o = optparse.OptionParser()
o.set_usage('omini_antplot_1chan.py [options] *omni_output.npz')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
a.scripting.add_standard_options(o,chan=True)
opts,args = o.parse_args(sys.argv[1:])
chan = int(opts.chan)

### Save Options ###
for filename in args:
    print filename,'-->',
    outfile = filename[:-4]+'.slice.npz'
    print outfile
    if os.path.exists(outfile): 
        print "     ... exists. skipping"
        continue
   
    meta,gains,vismdl,xtalk = from_npz(filename)
    slice_gains = {}
    slice_vismdl = {}
    slice_xtalk = {}
    for pol in gains:
        slice_gains[pol] = {}
        for ant,gain in gains[pol].iteritems():
            slice_gains[pol][ant] = n.array([g[chan] for g in gain])
    for pol in vismdl:
        slice_vismdl[pol] = {}
        for bl,vis in vismdl[pol].iteritems():
            slice_vismdl[pol][bl] = n.array([v[chan] for v in vis])
        slice_xtalk[pol] = {}
        for bl,xt in xtalk[pol].iteritems():
            slice_xtalk[pol][bl] = n.array(xt[chan])
    to_npz(outfile,meta,slice_gains,slice_vismdl,slice_xtalk)

