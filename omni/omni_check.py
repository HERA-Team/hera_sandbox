#! /usr/bin/env python

import omnical
import aipy
import pylab
import numpy as np
import capo
import pickle
import matplotlib.pyplot as plt
import optparse
import os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_check.py [options] *.npz')
aipy.scripting.add_standard_options(o,pol=True,max=True,drng=True)
o.add_option('--chisq',dest='chisq',default=False,action="store_true",
            help='Plot chisq.')
o.add_option('--gains',dest='gains',default=False,action="store_true",
            help='Plot gains of each antenna solved for.')
o.add_option('--chisqant',dest='chisqant',default=False,action="store_true",
            help='Plot chisqs per antenna.')
o.add_option('-i','--interactive',default=False,action="store_true",help='Launch IPython session before plotting stage')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',help='Path and name of calfile.')
o.add_option('--degen',dest='degen',action='store_true',help='Plot remaining array degeneracies. MUST also have the "gains" option on.')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


### Plot ChiSq ####
#if opts.pol == -1: pol = args[0].split('.')[3] #XXX hard-coded for *pol.npz files
pol = opts.pol
if opts.chisq == True:
    chisqs = []
    for i,file in enumerate(args):
        print 'Reading',file
        file = np.load(file)
        try: #reads *pol.npz files
            chisq = file['chisq '+str(pol)] #shape is (#times, #freqs)
        except: #reads .npz files
            chisq = file['chisq']
        for t in range(len(chisq)):
            chisqs.append(chisq[t])
            #chisq is sum of square of (model-data) on all visibilities per time/freq snapshot

    cs = np.array(chisqs)
    plt.imshow(np.log(cs),aspect='auto',interpolation='nearest',vmax=7,vmin=-6)
    plt.xlabel('Freq Channel',fontsize=10)
    plt.ylabel('Time',fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=8)
    plt.title('Omnical ChiSquare',fontsize=12)
    plt.colorbar()
    plt.show()


### Plot Gains or Chisqants ###
if opts.gains == True or opts.chisqant == True:
    gains = {} #or chisqant values, depending on option
    for i, f in enumerate(args): #loop over files
        print 'Reading',f
        file = np.load(f)
        for a in range(128):
            if opts.chisqant == True: #chisqant / gain
                try: value = file['chisq'+str(a)+pol[0]]/file[str(a)+pol[0]] #XXX only 0th element of pol
                except: continue
                try: gains[a].append(value)
                except: gains[a] = [value]
                vmax=0.05
                vmin=0.0
            if opts.gains == True:
                print "option doesn't exist yet"
                sys.exit()
                try: value = file[str(a)+pol[0]] #XXX only the 0th element of pol
                except: continue
                try: gains[a].append(value)
                except: gains[a] = [value]
                vmax=1.5
        file.close()
    for key in gains.keys():
        #gains[key] = np.vstack(numpy.abs(gains[key]))
        gains[key] = np.vstack(gains[key])
        mk = np.ma.masked_where(np.abs(gains[key]) == 1,np.abs(gains[key])).mask #flags
        gains[key] = np.ma.masked_array(gains[key],mask=mk) #masked array
    if opts.degen:
        #calculate degeneracies following Liu+'10 Eqns 11,12a,12b
        aa = aipy.cal.get_aa(opts.cal,np.array([0.15]))
        R = np.array([aa.get_baseline(0,i) for i in np.arange(112)]) # get positions relative to antenna 0
        Dx,Dy,Dz,Do = np.zeros_like(gains['0']),np.zeros_like(gains['0']),np.zeros_like(gains['0']),np.zeros_like(gains['0'])
        for i in np.arange(112):
            try:
                """
                I'm not assuming the array is planar, as they do in the paper. This is the "simple" version of
                non-coplanar degeneracy. See Liu+'10 Eqn 44 for the "real" version
                """
                Dx+=R[i,0]*np.angle(gains[str(i)])
                Dy+=R[i,1]*np.angle(gains[str(i)])
                Dz+=R[i,2]*np.angle(gains[str(i)])
                Do+=np.angle(gains[str(i)])
            except KeyError: continue

        if opts.interactive: import IPython;IPython.embed()

    #plotting stage
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if opts.degen:
        f,axarr = plt.subplots(1, 4,sharex=True,sharey=True)
        labels = ['x','y','z','overall']
        Dpars = [Dx,Dy,Dz,Do]
        if not opts.max is None: mx = opts.max
        else: mx = np.log10(np.abs(Dx)).max()
        if not opts.drng is None: mn = mx - opts.drng
        else: mn = np.log10(np.abs(Do)).min()
        for i in range(4):
            im = axarr[i].imshow(np.log10(np.abs(Dpars[i])),vmax=mx,vmin=mn,aspect='auto',interpolation='None')
            axarr[i].set_xlabel('Frequency bin')
            divider=make_axes_locatable(axarr[i])
            cax = divider.append_axes("right",size="20%",pad=0.05)
            cbar = plt.colorbar(im,cax=cax)
            axarr[i].set_title(labels[i]+' phase redundancy')
        axarr[0].set_ylabel('Integration number')
        #plt.tight_layout()
        #plt.subplots_adjust(top=0.85)    
        #plt.setp([a.get_yticklabels() for a in f.axes[1:]], visible=False)
    plt.show()
    
    #plot individual antenna gains
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
        plt.imshow(np.abs(gains[ant]),vmax=vmax,aspect='auto',interpolation='nearest')
        plt.title(ant,fontsize=10)
        plt.tick_params(axis='both',which='major',labelsize=6)
        plt.tight_layout()
        subplotnum += 1
    plt.show()
    print '1 sigma cut on median chisq: ',baddies
