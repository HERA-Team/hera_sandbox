#! /usr/bin/env python

import sys,os,optparse
import datetime
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pylab import *

o = optparse.OptionParser()
o.add_option('-o','--outfile',dest='outfile',default='today.png',help='Path to the output image file')
opts,args = o.parse_args(sys.argv[1:])

times,freqs,vTime,vFreq = None,None,[],None
flag_arr = []
for file in args:
    dat = np.load(file)
    files = dat.files
    dat['times']
    if times is None: times = dat['times']
    else: times = np.concatenate((times,dat['times']),axis=0)
    for i in range(len(dat['times'])):
        vTime.append(dat[str(i)].sum()/float(len(dat[str(i)])))
        if freqs is None: freqs = np.linspace(0.1,0.2,dat[str(i)].shape[0])
        if vFreq is None: vFreq = dat[str(i)].astype(np.float)
        else: vFreq += dat[str(i)]
        flag_arr.append(dat[str(i)])
vFreq /= vFreq.max()
flag_arr = np.array(flag_arr)
JD0 = np.floor(times[0])
JD = times - JD0

fig1 = figure(figsize=(8,10))
#plot v frequency
axf = fig1.add_subplot(311)
axf.plot(freqs,vFreq)
axf.set_xlabel('Frequency [GHz]')
#plot v time
axt = fig1.add_subplot(312)
axt.plot(JD,vTime,'.')
axt.set_xlabel('Days since JD %d'%JD0)
axwf = fig1.add_subplot(313)
axwf.imshow(flag_arr,aspect='auto',cmap='binary',extent=(freqs.min()*1e3,freqs.max()*1e3,times.max()*24,times.min()*24))
axwf.set_ylabel('time [hrs]')
axwf.set_xlabel('freq [MHz]')
#title
fig1.suptitle('RFI summary for JD'%(JD))
print 'saving to %s' % opts.outfile
fig1.savefig(opts.outfile,fmt='png')
