#! /usr/bin/env python

###
# NAME:    check_omnical.py
# AUTHOR:  Carina Cheng
# PURPOSE: Looks at Omnical outputs
###

import numpy
import matplotlib.pyplot as plt
import aipy
import omnical
import os, sys
import optparse


### Options ###
o = optparse.OptionParser()
o.set_usage('check_omnical.py [options]')
o.set_description(__doc__)
o.add_option('-f', '--freqchan', default=0,help='Frequency channel to plot. Default is 0.')
o.add_option('-t','--time',default=0,help='Time channel to plot. Default is 0.')
opts,args = o.parse_args(sys.argv[1:])


### Omnigain ###
if args[0][-8:] == 'omnigain':
    jds = []
    gains_per_freq = {}
    for omnifile in range(len(args)):
        print 'Reading',args[omnifile]
        dt = numpy.dtype([('jd',numpy.float64), ('ant', numpy.float32), ('nfreq', numpy.float32), ('gain',numpy.complex64,(203,))])
        file = numpy.fromfile(args[omnifile],dtype=dt)

        for i in range(len(file['jd'])):
            jd = float(file['jd'][i])
            ant = int(file['ant'][i])
            gains = file['gain'][i]
            gains_per_freq[jd,ant] = gains
        ants = numpy.unique(file['ant'])
        jds.append(numpy.unique(file['jd']))

    jds = numpy.array(numpy.concatenate(jds))
    subplot_count = 1
    plot_num = 1
    #print ' Bad Antennas:' 
    for a in range(len(ants)): #loop over antennas to plot
        #gain vs .time
        gs = numpy.abs([gains_per_freq[jd,ants[a]][:] for jd in jds])
            #gt is a 2D array with shape (#JDs, #freqs)
            #for an antenna, gt doesn't change along the time axis (but it does along the frequency axis
        #waterfall plot
        if subplot_count == 26:
            plot_num += 1
            subplot_count = 1
        plt.figure(plot_num,figsize=(10,10))
        plt.subplot(5,5,subplot_count)
        ymax=jds.max()-int(str(jds.min()).split('.')[0])
        ymin=jds.min()-int(str(jds.min()).split('.')[0])
        plt.imshow(gs,extent=(0,202,ymax,ymin),aspect='auto',interpolation='nearest',vmax=6)
        plt.xlabel('Freq Channel',fontsize=8)
        plt.ylabel('Time [JD]',fontsize=8)
        plt.tick_params(axis='both',which='major',labelsize=6)
        plt.text(0.05,ymin,"+%i"%int(str(jds.min()).split('.')[0]),fontsize=8)
        plt.title(ants[a],fontsize=10)
        #plt.colorbar()
        plt.tight_layout()
        subplot_count += 1

    """
        #1D plots
        g_med_f = numpy.median(gs,axis=1) #median gain along freq axis per time
        g_std_f = numpy.std(gs,axis=1) #stdev gain along freq axis per time
        gt = gs[:,opts.freqchan]/g_med_f 
        #print outliers 
        mean_gt = numpy.mean(gt)
        if numpy.abs(mean_gt-numpy.mean(g_med_f)) > numpy.mean(g_std_f):
            print '   '+str(ants[a])
        #gain vs. freq
        gf = numpy.abs(gains_per_freq[jds[opts.time],ants[a]])
        plt.figure(numpy.ceil(len(ants)/25.+1))
        plt.subplot(1,2,1)
        plt.plot(jds-int(str(jds[0]).split('.')[0]),gt)
        plt.xlabel('Time [JD]',fontsize=10)
        plt.ylabel('Normalized Gain for Freq Channel '+str(opts.freqchan),fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.text(0.92,-0.07,"+%i"%int(str(jds[0]).split('.')[0]),fontsize=8,transform=plt.gca().transAxes)
        plt.subplot(1,2,2)
        plt.plot(gf,label=ants[a])
        plt.xlabel('Freq Channel',fontsize=10)
        plt.ylabel('Gain',fontsize=10)
        plt.xlim(0,203)
        plt.tick_params(axis='both', which='major', labelsize=6)
    plt.legend(loc=1,prop={'size':7})
    plt.tight_layout()
    """  
    plt.show()
        

        
else:
    print 'Code not set up to read this file.'


