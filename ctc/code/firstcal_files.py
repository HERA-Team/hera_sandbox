#!/usr/bin/env python

import numpy 
import capo
import matplotlib.pyplot as plt
import os,sys
import optparse
from mpldatacursor import datacursor

o = optparse.OptionParser()
o.add_option('--gain',action="store_true",default=False, help='Plot the firstcal gains for every day for an antenna opts.ant')
o.add_option('-a','--ant',default=0, help='Antenna for which to plot gain solutions.')
o.add_option('--delay',action="store_true",default=False, help='Plot the firstcal delays vs. time for every antenna')
opts,args = o.parse_args(sys.argv[1:])

if opts.gain == False and opts.delay == False:
    print 'Specify --gain or --delay for plot.'
    exit()

npzs= args
pol = npzs[0].split('/')[-1].split('.')[3]
gains_ant = {}
gains_days = {}
delays_ant = {}
jds = []
for i,npz in enumerate(npzs):
    if opts.gain: m,g,v,x = capo.omni.from_npz(npz)
    if opts.delay: f = numpy.load(npz)
    jds.append('.'.join(npz.split('/')[-1].split('.')[1:3]))
    print 'Reading %s'%npz     
    for ant in range(112):
        #if ant == 16 or ant == 57 or ant == 24 : continue
        #if ant == 16 or ant == 95 or ant == 99 or ant == 100 or ant == 24 or ant == 29 or ant == 57: continue
        if opts.gain:
            try: gain = g[pol[0]][ant] #find gain
            except: gain = numpy.zeros_like(g[pol[0]][int(opts.ant)]) #make solutions 0 if doesn't exist (bad antenna)
            try: gains_ant[ant] = numpy.vstack((gains_ant[ant],gain)) #fill gains_ant
            except: gains_ant[ant] = gain
            day = int(npz.split('/')[-1].split('.')[1])
            try: gains_days[day][ant] = numpy.vstack((gains_days[day][ant],gain)) #fill gains_days
            except: 
                try: gains_days[day][ant] = gain
                except:
                    gains_days[day] = {}
                    gains_days[day][ant] = gain
        if opts.delay:
            try: delay = f['d'+str(ant)][0]
            except: delay = numpy.zeros_like(f['d'+str(opts.ant)][0])
            try: delays_ant[ant] = numpy.concatenate((delays_ant[ant],delay))
            except: delays_ant[ant] = delay
if opts.gain: #gain plot per day for one antenna
    f,axarr = plt.subplots(10,10,sharex=True,sharey=True)
    axs = axarr.ravel()
    for i,ax in enumerate(axs):
        try: 
            ax.tick_params(axis='both',which='both',labelsize=8)
            ax.set_title(gains_days.keys()[i],fontsize=10)
            toplot = gains_days[gains_days.keys()[i]][int(opts.ant)]
        except: continue
        ax.imshow(numpy.angle(toplot),aspect='auto')
    f.suptitle('Antenna '+str(opts.ant))#+' across Season 1')
    f.subplots_adjust(hspace=0.3)
    plt.show()
if opts.delay: #delay plot vs. time for all antennas
    for ant in delays_ant:
        plt.plot(delays_ant[ant],label=ant)
    datacursor(formatter='{label}'.format)
    plt.legend(loc=4,prop={'size':6})
    plt.xlabel('Time')
    plt.ylabel('Delay')
    plt.show()

"""
# Subplot per antenna 
f,axarr = plt.subplots(8,14,sharex=True,sharey=True)
axs = axarr.ravel()
for i,ax in enumerate(axs):
    ax.imshow(numpy.real(gains[i]),aspect='auto')
    ax.set_title(str(i)+pol[0],fontsize=10)
    ax.tick_params(axis='both',which='both',labelsize=8)
f.subplots_adjust(hspace=0.3)
plt.show()
"""    
