#! /usr/bin/env python

import sys,os,optparse
import datetime
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import dates
from pylab import *

o = optparse.OptionParser()
o.add_option('-o','--outfile',dest='outfile',default='today.png',help='Path to the output image file')
opts,args = o.parse_args(sys.argv[1:])

def fmt_times(jd):
    #get the y,m,d out of the jd
    ymd = int(np.floor(jd))
    #from aa.usno.navy.mil
    l = ymd+68569
    n = 4*l/146097
    l = l - (146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l - 1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l
    Y,M,D = i,j,k

    #get h,m,s out of the jd.
    jd -= np.floor(jd)
    jd = 24.*jd
    h = np.floor(jd)
    jd = 60.*(jd-h)
    m = np.floor(jd)
    s = np.floor(60.*(jd-m))

    #format it
    ifmt = '%Y/%m/%d %H:%M:%S'
    tstr = '%d/%d/%d %d:%d:%d'%(Y,M,D,h,m,s)
    t = time.strptime(tstr,ifmt)
    return datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec)

times,freqs,vTime,vFreq = None,None,[],None

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
vFreq /= vFreq.max()

ftimes = []
for t in times: ftimes.append(fmt_times(t)) 

fig1 = figure()
#plot v frequency
axf = fig1.add_subplot(211)
axf.plot(freqs,vFreq)
axf.set_xlabel('Frequency [GHz]')
#plot v time
axt = fig1.add_subplot(212)
axt.plot_date(times,vTime,'.')
#format dates
Mfmt = dates.DateFormatter('%m - %d')
mfmt = dates.DateFormatter('%H:%M')
axt.xaxis.set_major_locator(dates.DayLocator())
axt.xaxis.set_minor_locator(dates.HourLocator(interval=3))
axt.xaxis.set_minor_formatter(mfmt)
axt.xaxis.set_major_formatter(Mfmt)
axt.set_xlabel('UTC')
#title
f0 = ftimes[0]
fig1.suptitle('RFI summary for (Y-M-D) %d-%d-%d'%(f0.year,f0.month,f0.day))
print 'saving to %s' % opts.outfile
fig1.savefig(opts.outfile,fmt='png')
