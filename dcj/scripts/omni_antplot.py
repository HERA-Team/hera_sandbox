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

### Options ###
o = optparse.OptionParser()
o.set_usage('omini_antplot.py [options] *omni_output.npz')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
o.add_option('--maxsubplots',type=float,default=36,
        help='max number of plots to put on the same page [default 36')
a.scripting.add_standard_options(o,chan=True)
o.add_option('--nsigma',type=float,default=100,
    help='flag solutions with log(gain)>std(log(gain))*nsigma [default=100]')
o.add_option('--savemask',action='store_true')
o.add_option('--plot',action='store_true')
opts,args = o.parse_args(sys.argv[1:])

### Functions ###
def dB(x):
    return  10*n.ma.log10(n.ma.abs(x))

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future
rawgains = {}
gaintimes = {}
for filename in args:
    jd = file2jd(filename)
    print filename,jd
    meta,gains,vismdl,xtalk = from_npz(filename)
#    D = n.load(filename)
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
    gains[thing] = n.array(n.vstack(sortedgains))
#ants = [l.split(',')[2] for l in gains if l.find('gain')>-1]
print gains.keys()
antgains = gains.keys()
#antgains = [l for l in gains if l.find('gain')>-1]
print "found ",len(antgains),"antennas"
#xtalks = [l for l in gains if l.find('xtalk')>-1]

#do some flagging
for antgain in antgains:
    gains[antgain] = n.ma.masked_where(dB(gains[antgain])<-5,gains[antgain])
#rms over all data for each channel
antgainstd ={}
print "masking summary by antenna"
for antgain in antgains:
    antgainstd[antgain] = n.ma.std(dB(gains[antgain]))
    gains[antgain] = n.ma.masked_where(
                n.abs(dB(gains[antgain])-n.mean(dB(gains[antgain])))>antgainstd[antgain]*opts.nsigma,
        gains[antgain])


#use the median rms over all ants to do some stricter flagging
medianstd = n.median([antgainstd[antgain] for antgain in antgains])
print "median standard deviation",medianstd
for antgain in antgains:
    gains[antgain] = n.ma.masked_where(
        n.abs(dB(gains[antgain])-n.mean(dB(gains[antgain])))>(medianstd*opts.nsigma),gains[antgain])

if opts.savemask:
    print "save a mask file"
    for i in xrange(len(args)):
        masksave = {}
        for ant in n.sort(antgains):
            m = n.copy(gains[ant].mask)
            m.shape = n.array(rawgains[ant]).shape
            antpol = ant.split(',')[2]+ant.split(',')[0]
            masksave[antpol] = m[i].astype(float)
        maskfile = args[i].replace('.npz','.mask.npz')
        print maskfile
        n.savez(maskfile,**masksave)

####PLOT TIME###
if not opts.plot: sys.exit()
if not opts.chan is None:
    figure()
    for i,antgain in enumerate(antgains):
        plot(dB(gains[antgain])[:,opts.chan],label=antgain)
ylim([-10,10])
#plot things
ngains = len(antgains)
maxsubplots = min((opts.maxsubplots,ngains))
rows = n.ceil(n.sqrt(maxsubplots))
cols = n.ceil(maxsubplots/rows)
if True:
    plotnum=1
    for i,antgain in enumerate(antgains):
        if plotnum>maxsubplots:
            plotnum=1
            figure()
        subplot(rows,cols,plotnum)
        imshow(n.angle(gains[antgain]),interpolation='nearest',aspect='auto',vmin=0,vmax=2*n.pi)
        title(antgain)
        plotnum+=1
    plotnum=1
    figure()
    for i,antgain in enumerate(antgains):
        if plotnum>maxsubplots:
            plotnum=1
            figure()
        subplot(rows,cols,plotnum)
        imshow(dB(gains[antgain]),interpolation='nearest',aspect='auto',vmax=1,vmin=-1)
        title(antgain)
        plotnum+=1
    colorbar()
plotnum=1


show()
sys.exit()
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
    plt.imshow(numpy.log(cs),aspect='auto',interpolation='nearest',vmax=7,vmin=-6)
    plt.xlabel('Freq Channel',fontsize=10)
    plt.ylabel('Time',fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=8)
    plt.title('Omnical ChiSquare',fontsize=12)
    plt.colorbar()
    plt.show()


