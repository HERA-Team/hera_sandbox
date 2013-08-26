#! /usr/bin/env python
"""
Plot the amplitude, phase, and passband in an aipy cal file
Also a nice way to test that a file is being imported from the location you want. 

plot_cal.py mycal mycal2 etc
"""
import aipy as a, numpy as n,math as m
import sys, optparse
from pylab import *
o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--freqs',default='100_200',
    help='frequency range in MHz [defaul 100_200]')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])

def dB(x):
    return 10*n.log10(x)

aas = []
#generate the frequency axis from the input
freqs = n.linspace(float(opts.freqs.split('_')[0]),float(opts.freqs.split('_')[1]))
for cal in args:
    print "loading cal",cal
    aas.append(a.cal.get_aa(cal,freqs/1e3))
figure()
for aa in aas:
    subplot(221)
    prms= aa.get_params({'*':['dly','amp']})
    ants = prms.keys()
    ants.sort(key=lambda x:float(x))
    delays = [prms[ant]['dly'] for ant in ants]
    amps = [prms[ant]['amp'] for ant in ants]

    plot(ants,delays)
    xlabel('ant')
    ylabel('delay [ns]')
    subplot(222)
    plot(ants,dB(amps))
    xlabel('ant')
    ylabel('gain [dB]')
    subplot(223)
    plot(freqs,dB(aa[0].passband()))
    xlabel('freq [MHz]')
    ylabel('gain [dB]')
show()




