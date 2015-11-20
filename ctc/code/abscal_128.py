#! /usr/bin/env python

import aipy
import numpy
import capo
import os,sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy
import pylab
import optparse
import glob

##### 
# Fits a bandpass to 128 data, using 64 data
# Outputs *uvcRREOC files
#####


### Options ###
o = optparse.OptionParser()
o.set_usage('abscal_128.py *uvcRREO')
o.set_description(__doc__)
o.add_option('--plot',dest='plot',default=False,action="store_true",
            help='Plot bandpass fits.')
o.add_option('--abscal',dest='abscal',default=False,action="store_true",
            help='Apply absolute calibration to files.')
opts,args = o.parse_args(sys.argv[1:])


### Read in data and find matching LSTs ###

aa = aipy.cal.get_aa('psa898_v003',0.001,0.1,203) #parameters don't matter... only used to find LSTs
data128 = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_test/lstbin_fg/even/*.uv')) #LST-binned, bl-avged, FG-containing 128 data
data64 = numpy.sort(glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/lstbin_fg_even/*uvAV')) #LST-binned, bl-avged, FG-containing, abscal 64 data
#data128 = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_test/lstbin/even/*uvV')) #LST-binned, bl-avged 128 data
#data64 = numpy.sort(glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/lstbin_even_noxtalk/*uvGV')) #LST-binned, bl-avged 64 data

t_128 = []
d_128 = []
t_64 = []
d_64 = []
for file in data128: 
    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='all',polstr='xx',return_lsts=True)
    for tt,time in enumerate(t):
        t_128.append(time)
        d_128.append(d[d.keys()[0]]['xx'][tt])
for file in data64:
    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='all',polstr='I',return_lsts=True)
    for tt,time in enumerate(t):
        t_64.append(time)
        d_64.append(d[d.keys()[0]]['I'][tt])
d_128 = numpy.abs(d_128)
d_64 = numpy.abs(d_64)

all_factors = []
for t,time in enumerate(t_128):
    index1 = numpy.argmin(numpy.abs(t_64-time))
    time1 = t_64[index1]
    if numpy.abs(t_64[index1+1]-time) < numpy.abs(t_64[index1-1]-time):
        index2 = index1+1
        time2 = t_64[index2]
    elif numpy.abs(t_64[index1+1]-time) > numpy.abs(t_64[index1-1]-time):
        index2 = index1-1
        time2 = t_64[index2]
    else:
        time2 = time1
    if time1 < time2:
        interp = interp1d([time1,time2],[d_64[index1],d_64[index2]],axis=0)
        d_64_interp = interp(time)
    elif time1 > time2:
        interp = interp1d([time2,time1],[d_64[index2],d_64[index1]],axis=0)
        d_64_interp = interp(time)
    else:
        d_64_interp = d_64[index1]
    factors = d_64_interp/d_128[t]
    all_factors.append(factors)
    #plt.plot(numpy.abs(numpy.fft.ifft(factors,axis=1))**2,'.')
    plt.plot(factors,'.')
spectrum = numpy.median(all_factors,axis=0)
if opts.plot:
    #plt.plot(numpy.median(all_factors,axis=0))
    plt.xlabel('Frequency Channel')
    plt.ylabel('PSA64 / PSA128')
    plt.title('Ratios for every PSA128 integration')
    plt.show()


### Fit Bandpass ###
freqs = numpy.linspace(0,202,203)
bandpass = numpy.polyval(numpy.polyfit(freqs[numpy.isfinite(spectrum)],spectrum[numpy.isfinite(spectrum)],10),freqs)
if opts.plot:
    plt.plot(bandpass,'r-',label='polyfit')
    plt.plot(spectrum,'k-',label='median over all time')
    plt.xlabel('Frequency Channel')
    plt.ylabel('PSA64 / PSA128')
    plt.legend()
    plt.show()
    plt.plot(spectrum/bandpass,'b-',label='data divided by fit')
    plt.xlabel('Frequency Channel')
    plt.ylabel('PSA64 / PSA128')
    plt.legend()
    plt.show()
#scipy.optimize.curve_fit(f,freqs,spectrum)



### Calibrate Files ###

if opts.abscal:
    print 'Absolute Calibrating...'
    #mfunc here


