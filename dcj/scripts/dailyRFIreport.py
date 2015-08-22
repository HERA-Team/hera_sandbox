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
axf.plot(freqs*1e3,vFreq)
axf.set_xlabel('Frequency [MHz]')
#plot v time
axt = fig1.add_subplot(312)
local_hours = (JD*24.+2+12)%24
RD = times-1721424.5
datetimes = matplotlib.dates.num2date(RD) #mpl uses Rata Die dates
axt.plot_date(datetimes,vTime,'.',tz='Africa/Windhoek',)
axt.set_xlabel('SAST [hours]')
#axt.xaxis.set_major_locator(matplotlib.dates.HourLocator)
axt.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H"))
#axt.set_minor_locator()
axwf = fig1.add_subplot(313)
axwf.imshow(np.flipud(flag_arr),aspect='auto',interpolation='nearest',cmap='binary',extent=(freqs.min()*1e3,freqs.max()*1e3,RD[0],RD[-1]))

axwf.yaxis_date(tz='Africa/Windhoek')
axwf.invert_yaxis()
axwf.yaxis.set_major_locator(matplotlib.dates.HourLocator(interval=1))
axwf.yaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H"))


axwf.set_ylabel('time [hrs]')
axwf.set_xlabel('freq [MHz]')
#title
fig1.suptitle('RFI summary for %s (JD %d)'%(datetimes[0].strftime('%d-%b-%Y'),JD0))
if opts.outfile == 'today.png':
    opts.outfile = str(JD0)+'_RFI_summary_'+datetimes[0].strftime('%d-%b-%Y')+'.png'
print 'saving to %s' % opts.outfile
fig1.savefig(opts.outfile,fmt='png')
