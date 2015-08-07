#! /usr/bin/env python

###
# NAME:    check_omnical.py
# AUTHOR:  Carina Cheng
# PURPOSE: Looks at Omnical outputs (*omnigain, *omnichisq, *diagtxt)
# EXAMPLE: check_omnical.py *omnigain*
    # NOTE: Can only look at one type of file at one time
###

import numpy
import matplotlib.pyplot as plt
import aipy
import omnical
import capo
import pickle
import glob
import os, sys
import optparse


### Options ###
o = optparse.OptionParser()
o.set_usage('check_omnical.py [options]')
o.set_description(__doc__)
o.add_option('-f', '--freqchan', default=0,help='Frequency channel to plot (for 1D gain plot). Default is 0.')
o.add_option('-t','--time',default=0,help='Time channel to plot (for 1D gain plot). Default is 0.')
o.add_option('--mode',dest='mode',default='log',type='string',
            help='Mode to plot waterfall data in. Default is "log".')
opts,args = o.parse_args(sys.argv[1:])


### Omnigain ### 
#       -plots gains of each antenna solved for by Omnical (waterfall plots of time vs. frequency for each antenna)
#       -plots gains of each antenna vs. frequency and vs. time (2 plots)
if args[0][-8:] == 'omnigain':
    jds = []
    gains_per_freq = {}
    for omnifile in range(len(args)):
        print 'Reading',args[omnifile]
        dt = numpy.dtype([('jd',numpy.float64), ('ant', numpy.float32), ('nfreq', numpy.float32), ('gain',numpy.complex64,(203,))])
        file = numpy.fromfile(args[omnifile],dtype=dt)
        #jd is 1D array containing all times
        #ant is 1D array containing all antenna numbers, repeated as many times as there are 
        #gain is 2D array of shape (#times*#antennas, #freqs)
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
    print ' Bad Antennas: (>1 stdev away from median of all antennas)'
    all_gs = [] 
    for a in range(len(ants)): #loop over antennas to plot
        #gain vs .time
        gs = numpy.abs([gains_per_freq[jd,ants[a]][:] for jd in jds])
            #gt is a 2D array with shape (#JDs, #freqs)
            #for an antenna, gt doesn't change along the time axis (but it does along the frequency axis
        #waterfall plot
        all_gs.append(gs)
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

        #1D plots
        g_med_f = numpy.median(gs,axis=1) #median gain along freq axis per time
        g_std_f = numpy.std(gs,axis=1) #stdev gain along freq axis per time
        gt = gs[:,opts.freqchan]/g_med_f 
        #print outliers 
        mean_gt = numpy.mean(gt)
        #if numpy.abs(mean_gt-numpy.mean(g_med_f)) > numpy.mean(g_std_f):
            #print '   '+str(ants[a])
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
    gs_med = numpy.median(all_gs,axis=0) #median for all antennas
    gs_med = numpy.mean(gs_med)
    gs_std = numpy.std(all_gs,axis=0) #std for all antennas
    gs_std = numpy.mean(gs_std)
    for g,gs in enumerate(all_gs):
        if numpy.abs(numpy.mean(gs)-gs_med) > gs_std:
            print '   ',ants[g]
    plt.legend(loc=1,prop={'size':7})
    plt.tight_layout() 
    plt.show()

       
### ChiSquared ###
#       -plots chi-squared of Omnical as a waterfall plot (time vs. frequency)
elif args[0][-9:] == 'omnichisq':
    jds = []
    chisq_per_time = {}
    for omnifile in range(len(args)):
        print 'Reading',args[omnifile]
        dt = numpy.dtype([('jd',numpy.float64), ('nfreq', numpy.float32), ('chisq',numpy.float32,(203,))])
        file = numpy.fromfile(args[omnifile],dtype=dt)
        #jd is 1D array containing all times
        #nfreq is 1D array containing '203' repeated as many times as there are
        #chisq is 2D array of shape (# times, #freqs)
            #it changes with time and frequency
            #sum of square of (model-data) on all visibilities per time/freq snapshot
        for i in range(len(file['jd'])):
            jd = float(file['jd'][i])
            chisqs = file['chisq'][i]
            chisq_per_time[jd] = chisqs
            jds.append(jd)
    cs = [chisq_per_time[jd] for jd in jds]
    jds = numpy.array(jds)
    ymax = jds.max()-int(str(jds.min()).split('.')[0])
    ymin = jds.min()-int(str(jds.min()).split('.')[0])
    plt.imshow(numpy.log(cs),extent=(0,202,ymax,ymin),aspect='auto',interpolation='nearest',vmax=7,vmin=-6)
    plt.xlabel('Freq Channel',fontsize=10)
    plt.ylabel('Time [JD]',fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=8)
    plt.text(0.05,ymin,"+%i"%int(str(jds.min()).split('.')[0]),fontsize=10)
    plt.title('Omnical ChiSquare',fontsize=12)
    plt.colorbar()
    plt.show()
        

### Diagtxt Bad Antennas ###
#       -plots bad antenna counts as identified by Omnical
elif args[0][-7:] == 'diagtxt':
    ant_count = []
    for omnifile in range(len(args)):
        print 'Reading', args[omnifile]
        file = open(args[omnifile])
        for line in file:
            if line[0:7] == 'antenna':
                num = line[9:12]
                if num[2] == ',':
                    num = num[0:2]
                if num[1] == ',':
                    num = num[0]
                ant_count.append(int(num))
    n = plt.hist(ant_count,bins=128)
    n = n[0]
    ns = []
    for nn in range(len(n)):
        if n[nn] > len(args)/3: #1/3 of files
            ns.append(nn)
    print 'Bad antennas for at least 1/3 of the files:',ns
    plt.xlabel('Antenna Number')
    plt.ylabel('Number of files for which the antenna is bad')
    plt.title('Omnical Bad Antennas for '+str(len(args))+' Files')
    plt.show()


### First-Cal Applied Solutions from calpar File ###
#       -plots waterfall visibilities for selected baselines, corrected by the first_cal antenna solutions from the calpar.p file
elif args[0][:6] == 'calpar' and len(args) == 1:
    f = open(args[0],'rb')
    pol = args[0].split('.')[-2][-2:]
    file = numpy.array(pickle.load(f)[pol]) #shape = (#freqs, #ants)
    bls = ['55_27','60_39','42_37','33_6','52_7','21_15','41_19','53_21','67_47'] #baselines to look at
    compresseddata = glob.glob('/data4/paper/2014EoR/pot2/data2/2456926/zen.2456926.*.xx.uvcRRE')
    print 'Files =',compresseddata
    blstoget = ""
    for i in bls:
        blstoget += i+','
    blstoget = blstoget[:-1]
    print 'Reading files...'
    t,d,f = capo.arp.get_dict_of_uv_data(compresseddata,blstoget,polstr=pol)
    for b,bl in enumerate(bls):
        print 'Baseline '+str(b+1)+'/'+str(len(bls))
        ant1,ant2 = int(bl.split('_')[0]),int(bl.split('_')[1])
        calpar1,calpar2 = file[:,ant1],file[:,ant2]
        d_bl = d[aipy.miriad.ij2bl(int(bl.split('_')[0]),int(bl.split('_')[1]))][pol] #shape = (#times, #freqs)
        den = numpy.conj(calpar1)*calpar2
        d_cal = d_bl/den
        plt.subplot(3,3,b+1)
        capo.arp.waterfall(d_cal,mode=opts.mode)
        plt.colorbar()
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.ylabel('Time',fontsize=8)
        plt.xlabel('Freq Chan',fontsize=8)
        plt.title(bl,fontsize=10)
    plt.suptitle('first_cal-ed 128 File',fontsize=12,y=1.0)
    plt.tight_layout()
    plt.show()

     
else:
    print 'Code not set up to read this file.'


