#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import optparse, sys, os
import capo.uCal as uc

o = optparse.OptionParser()
o.set_usage('plot_uCal_result.py [.uCalResults.npz file]')
opts,args = o.parse_args(sys.argv[1:])
if len(args) == 0:
    print "WARNING!!! No arguments provided. Using testing parameters."
    resultsFile = './Data/zen.2456679.49577.xx.uCalResults.npz'
else: resultsFile = args[0]

dataFiles, allBandpasses, meanBandpass, stdBandpass, unflaggedChans, channelRMSs, overallChannelRMS, bandpassFit = uc.from_npz(resultsFile)

allUnflaggedChans = sorted(np.unique([item for sublist in unflaggedChans for item in sublist]))

#%% Plot bandpass
plt.figure(1); plt.clf()
for chans, bandpass in zip(unflaggedChans, allBandpasses):
    plt.plot(chans, np.abs(bandpass)[chans],':')
plt.errorbar(allUnflaggedChans,np.abs(meanBandpass)[allUnflaggedChans],yerr=np.abs(stdBandpass)[allUnflaggedChans], fmt='k')
plt.xlabel('channel'); plt.ylabel('abs(Bandpass)')

#%% Plot phase
plt.figure(2); plt.clf()
for chans, bandpass in zip(unflaggedChans, allBandpasses):
    plt.plot(chans, np.angle(bandpass)[chans],':')
plt.plot(allUnflaggedChans, np.angle(meanBandpass)[allUnflaggedChans],'k')
plt.xlabel('channel'); plt.ylabel('Bandpass Phase')

#%% Plot Channel RMS for each integration and all integrations
plt.figure(3); plt.clf()
for chans, channelRMS in zip(unflaggedChans, channelRMSs):
    plt.plot(chans, channelRMS[chans],':')
plt.plot(allUnflaggedChans, overallChannelRMS[allUnflaggedChans], 'k')
plt.xlabel('channel'); plt.ylabel('ChiSq/DoF')

plt.show()