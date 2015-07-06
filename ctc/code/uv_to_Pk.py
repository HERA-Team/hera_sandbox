#!/usr/bin/env python

"""

NAME:
     uv_to_Pk.py
PURPOSE:
     Plots P(k) from a UV file
EXAMPLE CALL:
     ./uv_to_Pk.py --file <path>
AUTHOR:
     Carina Cheng

"""

import numpy
import optparse
import capo
import aipy
import matplotlib.pyplot as plt
import os, sys

o = optparse.OptionParser()
o.set_usage('uv_to_Pk.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--file',dest='file', default='/home/cacheng/capo/ctc/tables/203files/pspec_Jy.uv', type='string',
             help='Path and name of UV file to plot.')
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is .1/203')
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
o.add_option('-c', dest='chan', default='95_115', type='string',
             help='Select channels to use. Example: 0_10 (includes channels 0 through 10). Default is 95_115.')
opts, args = o.parse_args(sys.argv[1:])

print 'Reading', opts.file

file = aipy.miriad.UV(opts.file)
data = []
times = []
for p,d,f in file.all(raw=True):
    pol = aipy.miriad.pol2str[file['pol']]
    if pol == 'xx':
        data.append(d) #length of d is number of freq channels
        times.append(p[1])

print 'Number of Time Integrations:', len(data) 
print 'Number of Frequency Channels:', len(data[0])

c1,c2 = opts.chan.split('_')
c1 = int(c1)
c2 = int(c2)

data = numpy.array(data)
times = numpy.array(times)
#print numpy.average(numpy.abs(data)**2,axis=0)
#factors

freqs = numpy.arange(opts.sfreq,opts.sfreq+(opts.nchan*opts.sdf),opts.sdf)
freqs_subset = freqs[c1:c2]
cent_freq = numpy.average(freqs_subset)

etas = capo.pspec.f2eta(freqs_subset)
z = capo.pspec.f2z(cent_freq)
ks = capo.pspec.dk_deta(z)*etas
 
X2Y = capo.pspec.X2Y(z)
B = freqs_subset[-1]-freqs_subset[0]
beam = ((numpy.polyval(capo.pspec.DEFAULT_BEAM_POLY,freqs_subset)))*2.35
factor = capo.pspec.jy2T(freqs_subset)

windowtype = 'hamming'
data = data[:,c1:c2] #numpy.transpose(numpy.transpose(data)[c1:c2])
data = data*factor
window = aipy.dsp.gen_window(data.shape[1],windowtype)
window.shape = (1,window.size)
data = data*window

#B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] # this is wrong if we aren't inverting
# the window post delay transform (or at least dividing out by the gain of the window)
# For windowed data, the FFT divides out by the full bandwidth, B, which is
# then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand

B = B*capo.pfb.NOISE_EQUIV_BW[windowtype] #need to change B because of windowing

ps_data = numpy.abs(numpy.fft.ifft(data,axis=1))**2 #(FT freq-axis)^2 
print 'Number of Channels Used:', len(ps_data[0])

"""
#bootstrapping

NBOOT=100

print 'Bootstrapping...'

final_ps = []
for boot in range(NBOOT):
    boot_ps = []
    if boot % 10 == 0: print boot
    for t in range(len(data)): #draws a number of samples = number of times
        t = numpy.random.choice(range(numpy.array(ps_data).shape[0]),replace=True)
        ps = ps_data[t]
        #print len(ps) #number of freqs               
        boot_ps.append(ps)
    #print len(boot_ps) #number of times
    median_ps = numpy.median(boot_ps,axis=0)
    final_ps.append(median_ps)

print 'Sorting bootstraps...'
"""

final_ps = [numpy.average(ps_data,axis=0)] #simplest case (no bootstrapping)

#print len(final_ps) #number of boots
#print len(final_ps[0]) #number of freqs in first boot sample

pk = numpy.average(final_ps,axis=0) #average Pk per k
scalar = X2Y*B*beam
pk = pk*scalar
#pk_boot = numpy.sort(numpy.array(final_ps).real,axis=0)
err = numpy.std(pk,axis=0) 

#plt.errorbar(ks,pk,yerr=err,fmt='k.')
#plt.show()

#folding k's

folded_pk = []
folded_ks = []
for i in range(int(len(ks)/2.)+1):
    if i==0:
        folded_pk.append(pk[i])
        folded_ks.append(ks[i])
    elif len(ks)%2==0 and i==len(ks)/2:
        folded_pk.append(pk[i])
        folded_ks.append(numpy.abs(ks[i]))
    else:
        low_ind = i
        high_ind = len(ks)-i
        folded_pk.append((pk[low_ind]+pk[high_ind])/2.)
        folded_ks.append(ks[low_ind])

#plotting

plt.subplot(121)
plt.gca().set_yscale('log')
plt.plot(ks,pk,'k.')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.grid()
#plt.ylim(1e5,3e16)

plt.subplot(122)
plt.gca().set_yscale('log')
plt.plot(folded_ks,((numpy.array(folded_ks)**3)/(2*numpy.pi**2))*folded_pk,'k.')
plt.xlabel('k')
plt.ylabel('$\Delta^2$(k)')
plt.grid()
plt.ylim(1e0,1e9)

plt.show()
