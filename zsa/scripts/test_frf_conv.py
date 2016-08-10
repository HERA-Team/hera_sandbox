#!/usr/bin/env python
import aipy as a, numpy as n
import capo, optparse, sys
import frf_conv
import pylab as p

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from\
    flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--xtalk', type='float', default=-1,
    help='Minimum fringe rate (in Hz) to allow.  Anything varying slower than\
    this is considered crosstalk.  A negative value indicates nothing should be\
    considered crosstalk.  Default is -1.  A typical value for PAPER is 6e-5.')
o.add_option('--max_fr_frac', type='float', default=1.,
    help='Fractional width of the fringe-rate filter range to extract.Default\
    1.') 
o.add_option('--min_fr_frac', type='float', default=-.3,
    help='Fractional width of the fringe-rate filter range to extract.Default\
    -.3') 
o.add_option('--ntaps', type='int', default=20,
    help='Number of taps to use in FIR filter.')
o.set_usage('fringe_rate_filter.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
print "inttime = uv['inttime']*4 #hack for *E files:inttime set incorrectly"
inttime = uv['inttime']*4 #hack for *E files:inttime set incorrectly
nchans = uv['nchan']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

#kernels for all frequencies
#beam_w_fr= frf_conv.get_beam_w_fr(aa, (49,41))#, binwidth=5e-6)   
beam_w_fr= frf_conv.get_beam_w_fr(aa, (1,4))#, binwidth=5e-6)   
tb, ks = frf_conv.get_fringe_rate_kernels(beam_w_fr, 42.8, 14*14)
p.figure(100)

p.plot(tb,ks[160].real)
p.plot(tb,ks[160].imag)

#capo.arp.waterfall(beam_w_fr, mx=-1, drng=5)

#take = n.where(n.log10(beam_w_fr) > -5)
#maxtake = n.max(take[1]) 
#mintake = n.min(take[1])
#print maxtake, mintake

#bins = bins[mintake:maxtake]
#beam_w_fr = beam_w_fr[:,mintake:maxtake]

#kernels = n.fft.ifft(beam_w_fr)
#kernels = n.fft.fftshift(n.fft.ifft(beam_w_fr, axis=-1), axes=-1)

import IPython
IPython.embed()

files = args
#Get times, data, and flags. Number of integrations.
times, data, flags = capo.arp.get_dict_of_uv_data(files,'49_41','I',verbose=1)
nints = len(times); print nints
mir_bl = a.miriad.ij2bl(49,41)
data = data[mir_bl]['I']
print len(data[0])
#simulate data
data = n.random.randn(nints,len(data[0])) + 1j*n.random.randn(nints,len(data[0]))
data = data.transpose()
flags = n.zeros_like(data)


#window function in time domain.
window = a.dsp.gen_window(len(bins),'blackman-harris')
print len(window)
#p.figure(2)
s_rate = bins[2] - bins[1]
fft_times = n.fft.fftfreq(bins.size, s_rate)
fft_times = n.fft.fftshift(fft_times)

d = data
w = n.logical_not(flags).astype(float)
p.figure(1)
capo.arp.waterfall(beam_w_fr, mx=-5, drng=3); p.colorbar(shrink=.5)
p.title('Gassian beam filter (fr-rate)')


#get convolution kernel. Note kernels should be in time order.
kernels = kernels * window

p.figure(2)
capo.arp.waterfall(kernels, mx=-5, drng=3);p.colorbar(shrink=.5)
p.title('convolution kernel (time)')
p.figure(10)
capo.arp.waterfall(n.fft.ifft(kernels), mx=-5, drng=3) ; p.colorbar(shrink=.5)
p.title('convolution kernel (fr_rate)')

p.figure(3)
p.imshow(n.abs(d))
p.title('input data ( amp vs. time for single freq.)')
print n.mean(d)
print n.std(d)

p.figure(4)
_d = []
_w = []
for kf,df,wf in zip(kernels, d, w):
    _d.append( n.convolve(df, kf, mode='same') )
    _w.append( n.convolve(wf, n.abs(kf), mode='same') )

_d = n.array(_d)
_w = n.array(_w)
#p.plot(n.abs(kernel)/n.max(n.abs(kernel)))
filtered = _d/_w
p.figure(4); capo.arp.waterfall(_d); p.colorbar(shrink=.5)
p.title('convolved data')
p.figure(7); capo.arp.waterfall(_w); p.colorbar(shrink=.5)
p.title('convolved wgts')
p.figure(8); capo.arp.waterfall(filtered); p.colorbar(shrink=.5)
p.title('convolved (time)')
p.figure(5)
capo.arp.waterfall(n.fft.ifft(filtered)); p.colorbar(shrink=.5)
p.title('convolved (freq)')


p.show()
