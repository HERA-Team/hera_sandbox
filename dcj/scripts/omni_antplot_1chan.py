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
from capo.omni import from_npz
### Functions ###
def dB(x):
    return  10*n.ma.log10(n.ma.abs(x))
### Options ###
o = optparse.OptionParser()
o.set_usage('omni_antplot_1chan.py [options] *omni_output.npz')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
o.add_option('--reload',action='store_true',
    help='By default this script looks for a cache to replot. this option triggers a re-read of the data')
a.scripting.add_standard_options(o,chan=True)
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
if (not os.path.exists('omni_antplot_1chan-cache.npz')) or opts.reload:
    pol = 'xx' #XXX shouldn't be hard-coded in the future
    rawgains = {}
    gaintimes = {}
    for filename in args:
        jd = file2jd(filename)
        print filename,jd
        meta,gains,vismdl,xtalk = from_npz(filename)
        for pol in gains:
            for ant,gain in gains[pol].iteritems():
                antpol= str(ant)+pol
                try:
                    rawgains[antpol].append(gain)
                    gaintimes[antpol].append(jd)
                except(KeyError):
                    rawgains[antpol] = [gain]
                    gaintimes[antpol] = [jd]
    #sort the things by jd
    gains  ={}
    for thing in rawgains:
        jdorder = n.argsort(gaintimes[thing])
        sortedgains = [rawgains[thing][i] for i in jdorder]
        gains[thing] = n.concatenate(sortedgains)
    n.savez('omni_antplot_1chan-cache',**gains)
else:
    gains = n.load('omni_antplot_1chan-cache.npz')
antgains = gains.keys()
print "found ",len(antgains),"antennas"
from mpldatacursor import datacursor
for i,antgain in enumerate(antgains):
    plot(dB(gains[antgain]),label=antgain)
datacursor(formatter='{label}'.format)
#legend()
show()
