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
o.set_usage('omni_check.py [options] *.npz')
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--chisq',dest='chisq',default=False,action="store_true",
            help='Plot chisq.')
o.add_option('--gains',dest='gains',default=False,action="store_true",
            help='Plot gains of each antenna solved for.')
o.add_option('--chisqant',dest='chisqant',default=False,action="store_true",
            help='Plot chisqs per antenna.')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
opts,args = o.parse_args(sys.argv[1:])


### Plot ChiSq ####
if opts.chisq == True:
    if opts.pol == -1:
        pol = args[0].split('.')[3] #XXX hard-coded for *pol.npz files
    chisqs = []
    for i,file in enumerate(args):
        print 'Reading',file
        file = numpy.load(file)
        try: #reads *pol.npz files
            chisq = file['chisq '+str(pol)] #shape is (#times, #freqs)
        except: #reads .npz files
            chisq = file['chisq']
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
if opts.gains == True or opts.chisqant == True:
    gains = {} #or chisqant values, depending on option
    for i, file in enumerate(args): #loop over files
        print 'Reading',file
        file = numpy.load(file)
        for key in file.keys(): #loop over antennas
            gain = file[key]
            if key[0] != '<' and key[0] != '(' and key[0].isalpha() != True and opts.gains == True:
                antnum = key[:-1]
                try: gains[antnum].append(gain)
                except: gains[antnum] = [gain]
                vmax=1.5
            if key[0] == 'c' and opts.chisqant == True and len(key) > 5: #if plotting chisq per ant
                antnum = key.split('chisq')[1][:-1]
                try: gains[antnum].append(gain)
                except: gains[antnum] = [gain]
                vmax=0.5
    for key in gains.keys():
        gains[key] = numpy.vstack(numpy.abs(gains[key])) #cool thing to stack 2D arrays that only match in 1 dimension
        mk = numpy.ma.masked_where(gains[key] == 1,gains[key]).mask #flags
        gains[key] = numpy.ma.masked_array(gains[key],mask=mk) #masked array
    #Plotting
    f,axarr = plt.subplots(8,14,figsize=(14,8),sharex=True,sharey=True)
    for ant in range(max(map(int,gains.keys()))+1):
        i1 = ant/14 #row number
        i2 = ant%14 #col number
        axarr[i1,i2].set_title(ant,fontsize=10)
        axarr[i1,i2].tick_params(axis='both',which='both',labelsize=8)
        try: axarr[i1,i2].imshow(gains[str(ant)],vmax=vmax,aspect='auto',interpolation='nearest')
        except: continue
    f.subplots_adjust(hspace=0.7)
    plt.show()
