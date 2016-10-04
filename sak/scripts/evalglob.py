#!/usr/bin/env python

import aipy, sys, capo, numpy as np
from matplotlib import pylab

def get_bl_avgs(d,t0=0,t1=14):
    """
    Takes time average of all BLs individually
    assumes d contains only one type of BL
    """
    bl_avgs = np.zeros((len(d.keys()),203))
    for i,bl in enumerate(d.keys()):
        tavg = np.nanmean(d[bl]['xx'][t0:t1,:],axis=0) #time avg of single BL
        bl_avgs[i,:] = tavg
    return bl_avgs

pp = 'xx'
print 'Reading'
tinfo,d,f = capo.arp.get_dict_of_uv_data(sys.argv[1:],antstr='cross',polstr=pp)
print 'done'

T = d[d.keys()[0]][pp].shape[0]
fbins = range(203)
freqs = np.linspace(.1,.2,num=203)
delays = np.fft.fftfreq(203,d=np.diff(freqs)[0])
weights = aipy.dsp.gen_window(203,'blackman-harris')

final_wf, final_df = np.zeros_like(d[d.keys()[0]][pp]), np.zeros_like(d[d.keys()[0]][pp])
bad_bl = [(50,54),(56,59),(34,51),(0,26),(0,44),(30,50),(12,26),(0,7),(0,46),(12,50),(30,59),(0,38),(0,60)]
for bl in d.keys():
    if bl in bad_bl: continue 
    final_wf += d[bl][pp]

final_wf /= len(d.keys())

for t in range(final_wf.shape[0]): 
    _d,info = aipy.deconv.clean(np.fft.ifft(final_wf[t,:]*weights), np.fft.ifft(np.ones((203))),tol=1e-3)
    final_df[t,:] = np.fft.fftshift(np.abs(_d+info['res']))

pylab.imshow(np.log10(np.abs(final_wf)),aspect='auto',interpolation='nearest');pylab.show()
pylab.close()   
pylab.imshow(np.log10(np.abs(final_df)),aspect='auto',interpolation='nearest');pylab.show()
pylab.close()

av=np.nanmean(final_df[200:800,:],axis=0)
pylab.plot(np.fft.fftshift(delays),av)
pylab.plot(np.fft.fftshift(delays)[102:],0.5*(np.flipud(av[:101])+av[102:]),'g--',lw=2)
pylab.show()


#import IPython;IPython.embed()
