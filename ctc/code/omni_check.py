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
o.set_usage('omni_check.py [options] *omni_output.npz')
o.add_option('--chisq',dest='chisq',default=False,action="store_true",
            help='Plot chisq.')
o.add_option('--gains',dest='gains',default=False,action="store_true",
            help='Plot gains of each antenna solved for.')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future


### Plot ChiSq ####
if opts.chisq == True:
    chisqs = []
    for i,file in enumerate(args):
        print 'Reading',file
        file = numpy.load(file)
        chisq = file[str(pol)+',chisq'] #shape is (#times, #freqs)
        for t in range(len(chisq)):
            chisqs.append(chisq[t])
        #chisq is sum of square of (model-data) on all visibilities per time/freq snapshot
    cs = numpy.array(chisqs)
    plt.imshow(numpy.log(cs),aspect='auto',interpolation='nearest',vmax=7,vmin=-6)
    plt.xlabel('Freq Channel',fontsize=10)
    plt.ylabel('Time',fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=8)
    plt.title('Omnical ChiSquare',fontsize=12)
    plt.colorbar()
    plt.show()


### Plot Gains ###
if opts.gains == True:
    gains = {}
    for i, file in enumerate(args): #loop over files
        print 'Reading',file
        file = numpy.load(file)
        for key in file.keys(): #loop over antennas
            if "gains" in key:
                gain = file[key]
                antnum = key.split(',')[2]
                try: gains[antnum].append(gain) 
                except: gains[antnum] = [gain]
    for key in gains.keys():
        gains[key] = numpy.vstack(gains[key]) #cool thing to stack 2D arrays that only match in 1 dimension
        print numpy.array(gains[key]).shape

    subplotnum = 1
    plotnum = 1
    plt.figure(plotnum,figsize=(10,10))
    for ant in gains.keys(): #loop over antennas
        if subplotnum == 26:
            #break #only generate one page of plots (faster for testing) 
            plotnum += 1
            plt.figure(plotnum,figsize=(10,10))
            subplotnum = 1
        plt.subplot(5,5,subplotnum)
        plt.imshow(numpy.abs(gains[ant]),vmax=1.5,aspect='auto',interpolation='nearest')
        plt.title(ant,fontsize=10)
        plt.tick_params(axis='both',which='major',labelsize=6)
        subplotnum += 1
    plt.tight_layout()
    plt.show()

