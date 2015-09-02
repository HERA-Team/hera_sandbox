#! /usr/bin/env python

import omnical
import aipy
import pylab
import numpy
import capo
import pickle
import matplotlib.pyplot as plt
import optparse
import os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omini_check.py [options] *omni_output.npz')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future


### Plot ChiSq ####
if args[0][-3:] == 'npz':
    chisqs = []
    for i,file in enumerate(args):
        print 'Reading',file
        file = numpy.load(file)
        chisq = file[str(pol)+',chisq'] #shape is (#times, #freqs)
        for t in range(len(chisq)):
            chisqs.append(chisq[t])
        #chisq is sum of square of (model-data) on all visibilities per time/freq snapshot
    cs = numpy.array(chisqs)
    plt.imshow(numpy.log(cs),aspect='auto',interpolation='nearest')#,vmax=7,vmin=-6)
    plt.xlabel('Freq Channel',fontsize=10)
    plt.ylabel('Time',fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=8)
    plt.title('Omnical ChiSquare',fontsize=12)
    plt.colorbar()
    plt.show()


