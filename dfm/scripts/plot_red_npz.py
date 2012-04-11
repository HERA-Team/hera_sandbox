#! /usr/bin/env python 

import aipy as a
import numpy as np
import os,sys,optparse
from pylab import *

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

time,amp,phs = None,None,None

for npzfile in args:
    print 'Reading', npzfile
    F = np.load(npzfile)
    if time is None:
        time = F['JD']
        amp = F['GAIN']
        phs = F['PHS']
    else:
        time = np.concatenate((time,F['JD']))
        amp = np.concatenate((amp,F['GAIN']))
        phs = np.concatenate((phs,F['PHS']))

time -= float(int(time[0])/1000)*1000.

man_flag_chans = a.scripting.parse_chans('0_130,755_777,1540,1704,1827,1868,1901_2047',time.shape[0])
flags = np.zeros_like(amp)
flags[:,man_flag_chans] = 1
flags = np.where(np.abs(amp) >= 10.,1,flags) 

amp = np.ma.array(amp,mask=flags)
phs = np.ma.array(phs,mask=flags)

figure(0)

subplot(121)
imshow(amp,aspect='auto')
colorbar()

subplot(122)
imshow(phs,aspect='auto')
colorbar()

draw()

figure(1)

subplot(211)
plot(time,np.mean(amp,axis=1))
ylabel('Amplitude ratio [db]')
xticks([])

subplot(212)
plot(time,np.mean(phs,axis=1))
xlabel('Julian Day')
ylabel('Phase Difference [rad]')

subplots_adjust(hspace=0)

draw()

show()

