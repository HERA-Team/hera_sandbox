import aipy, sys, capo, numpy as np
from matplotlib import pylab

DELAY = True
if DELAY: offset=35
else: offset=1000

FINAL_AVG = False

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

pieces = np.power(2,np.array(range(6)))
weights = aipy.dsp.gen_window(203,'blackman-harris')

for i,num in enumerate(pieces):
    final_avgs = np.zeros((num,203))
    times = np.linspace(0,T,num=num+1,dtype='int')
    for j,t in enumerate(times):
        try:
            t0,t1 = t,times[j+1]
            bl_avgs = get_bl_avgs(d,t0=t0,t1=t1)
            final_avgs[j,:] = np.nanmean(bl_avgs,axis=0)
        except IndexError:
            break
    for k in range(final_avgs.shape[0]): 
        if not DELAY: pylab.plot(fbins,final_avgs[k,:]+(offset*k),'b-', alpha=0.2)
        else:
            _d,info = aipy.deconv.clean(np.fft.ifft(final_avgs[k,:]*weights), np.fft.ifft(np.ones((203))),tol=1e-3) #XXX is the ifft(ones) correct?
            pylab.plot(np.fft.fftshift(delays), np.fft.fftshift(np.abs(_d+info['res']))+(offset*k), 'b-', alpha=0.5)
    
        
    favg=np.nanmean(final_avgs,axis=0)
    fstd=np.nanstd(final_avgs,axis=0)
    if not DELAY:
        if FINAL_AVG:
            pylab.plot(fbins,favg,'k-',lw=2)
            pylab.fill_between(fbins,favg-fstd,favg+fstd,alpha=0.1)
        pylab.xlim(0,203)
        pylab.xlabel(r'Frequency bin',size=15)
    else:
        _d_a,info_a = aipy.deconv.clean(np.fft.ifft(favg*weights), np.fft.ifft(np.ones((203))),tol=1e-3) #XXX is the ifft(ones) correct?
        if FINAL_AVG:
            pylab.plot(np.fft.fftshift(delays), np.fft.fftshift(np.abs(_d_a+info_a['res'])), 'k-', lw=2)
        #pylab.ylim(0,15)
        pylab.xlabel(r'Delay (ns)',size=15)
        pylab.xlim(-500,500) 
        pylab.axvline(x=-100)
        pylab.axvline(x=100)
    pylab.ylabel(r'Re(Vis) (arb).',size=15)
    pylab.suptitle(r'xtalk$_{%i}$'%(num),size=15)
    pylab.show()
    
