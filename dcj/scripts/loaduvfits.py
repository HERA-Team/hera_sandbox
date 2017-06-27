from astropy.io import fits
import os,sys
import numpy as n
from pylab import *
import capo
import aipy as a
poles = ['xx','yy','xy','yx']
WINDOW='blackman-harris'
def i2a(i):
    #convert from MWA to unit ant index
    tens = int(i/8)+1
    ones = i%8
    return tens*10+ones
def a2i(a):
    #convert from unit to MWA ant index
    eights = int(a/10)-1
    ones = a%10
    return eights*8+ones
def gethdufreqs(D):
    nfreqs = D.header['NAXIS4']
    f0  = D.header['CRVAL4']
    df = D.header['CDELT4']
    fi0 = D.header['CRPIX4']
    return n.arange(f0-df*fi0,f0-df*fi0+nfreqs*df,df)
def length(A):
    return n.sqrt(n.dot(A,A))
def gentle(dd, _w, tol=1e-4, area=None, stop_if_div=False,
                    maxiter=100,verbose=False):
    inside_res = n.std(dd[area!=0])
    outside_res = n.std(dd[area==0])
    print inside_res,outside_res
    cc = n.zeros_like(dd)
    while(inside_res>outside_res):
        _d_cl, info = a.deconv.clean(dd, _w, tol=1e-1, area=area, stop_if_div=False,
                    maxiter=100,verbose=False)
        res = info['res']
        inside_res = n.std(res[area!=0])
        outside_res = n.std(res[area==0])
        dd = info['res']
        cc += _d_cl
        print inside_res,outside_res
    return cc,info
for uvfile in sys.argv[1:]:
    F = fits.open(uvfile)
    D = F[0]
    times = D.data.field('DATE')
    t_int = n.diff(times).max()*24*3600
    print "t_int = ",t_int
    totalt = n.ceil((times.max() - times.min())*24*3600 + t_int)
    print "total observing time = ",totalt
    ntimes = n.ceil(totalt/t_int)
    print "ntimes = ",ntimes
    bls = D.data.field('BASELINE')
    ant2 = (bls%256).astype(int)
    ant1 = ((bls-ant2)/256).astype(int)
    uvws  = n.array(zip(D.data.field('UU'),D.data.field('VV'),D.data.field('WW')))
    Nants = len(set(ant2))
    Nmax = ant2.max()
    freqs = gethdufreqs(D)
    df = n.diff(freqs)[0]
    print "n_ant = ",Nants
    print n.sum((ant2-ant1)==0)
    print "forming up i,j,t,chan,pol matrix"
    Nbl = D.data.field('DATA').shape[0]
    print Nbl
    Nfreqs = D.data.field('DATA').shape[3]
    Npol = D.data.field('DATA').shape[4]
    if False:
        #plot autos
        for i in range(Nbl):
            if ant1[i]!=ant2[i]:continue
            data = D.data.field('DATA')[i,:,:,:,:,0] + 1j*D.data.field('DATA')[i,:,:,:,:,0]
            data = data.squeeze()
            plot(n.abs(data[:,0]),label=str(ant1[i]))
    if False:
        gains = n.zeros((Nmax,Nfreqs,2))
        count = n.zeros((Nmax,Nfreqs,2))
        for i in range(Nbl):
            if ant1[i]==ant2[i]:
                data = D.data.field('DATA')[i,:,:,:,:2,0] + 1j*D.data.field('DATA')[i,:,:,:,:2,0] 
                mask = D.data.field
                data = data.squeeze()
                gains[ant1[i]-1] += data
                count[ant1[i]-1] += 1.
        gains = n.ma.masked_invalid(n.array(gains)/n.array(count))
        gains /= n.ma.mean(gains)
        antgains = n.ma.mean(gains,axis=1).squeeze()
        antgains_err= n.ma.std(gains,axis=1).squeeze()

        changains = n.ma.mean(gains,axis=0).squeeze()
        changains_err = n.ma.std(gains,axis=0).squeeze()

        subplot(211)
        print antgains[:,0].shape 
        errorbar(range(len(antgains)),antgains[:,0],yerr=antgains_err,label='xx')
        errorbar(range(len(antgains)),antgains[:,1],yerr=antgains_err,label='yy')
        subplot(212)
        errorbar(range(len(changains)),changains[:,0],yerr=changains,label='xx')
        errorbar(range(len(changains)),changains[:,1],yerr=changains,label='yy')
        show()
    if True:
        window = None
        for i in range(Nbl):
            if ant1[i]==a2i(11) and ant2[i]==a2i(12):
                data = D.data.field('DATA')[i,:,:,:,0,0] + 1j*D.data.field('DATA')[i,:,:,:,0,1] 
                mask = (D.data.field('DATA')[i,:,:,:,0,2]==0).squeeze()
                print n.sum(mask)/float(mask.size)
                _w = n.fft.ifft(1-mask)
                figure()
                fill_between(freqs/1e6,1-mask.squeeze())
                ylim([0,2])
                show()
                data = n.ma.masked_where(mask,data.squeeze())
                uvw = uvws[i]
                print "bl length"
                print length(uvw)
                print 'delay'
                bldelay = length(uvw)*1e9
                print bldelay
                print 'uvw=',uvw               
                sys.exit()
                if window is None:
                    window = a.dsp.gen_window(data.shape[0],WINDOW)
                dd      =   n.fft.ifft(window*data)
                delays  =   n.fft.fftfreq(data.shape[0],df/1e9) #delays in ns
                
                area = n.zeros_like(dd).astype(n.int)
                area[n.abs(delays)<bldelay] = 1
                print "area"
                print n.sum(area)
                print "nan count"
                print n.sum(n.isnan(dd)),n.sum(n.isnan(_w))
                gain = n.abs(_w[0])
                print "gain = ",gain
                #_w = n.abs(_w)
                if gain>0:
                    #_d_cl, info = a.deconv.clean(dd, _w, tol=1e-4, area=area, stop_if_div=False,
                    #maxiter=100,verbose=False)
                    _d_cl, info = gentle(dd, _w, tol=1e-4, area=area, stop_if_div=False,
                    maxiter=100,verbose=False)

                    #_d_cl += info['res']
                else:
                    _d_cl = dd
                    continue
                dd = n.fft.fftshift(dd)
                delays = n.fft.fftshift(delays)
#clean
                subplot(311)
                plot(delays,n.fft.fftshift(n.abs(_w)))
                subplot(312)


                plot(delays,n.abs(dd))
                vlines([-1*bldelay,bldelay],[n.abs(dd).min()]*2,[n.abs(dd).max()]*2)
                subplot(313)
                if gain>0:
                    _d_cl = n.fft.fftshift(_d_cl)
                    res = n.fft.fftshift(info['res'])
                    plot(delays,n.abs(dd),'k')
                    plot(delays,n.abs(_d_cl),'b')
                    plot(delays,n.abs(res),'r')
                    plot(delays,n.abs(_d_cl +res),'g')
                    vlines([-1*bldelay,bldelay],[n.abs(dd).min()]*2,[n.abs(dd).max()]*2)
                    xlim([-2*bldelay,2*bldelay])
                show()
                #clean!!!
#                area = n.ones(dd.size, dtype=n.int)
#                area[uthresh:lthresh] = 0
#                _d_cl, info = a.deconv.clean(dd, _w, tol=opts.clean, area=area, stop_if_div=False, maxiter=100)
#        for i in list(set(ant1)):
#            plot(antgains[i],label=str(i))
    legend()
    show()
