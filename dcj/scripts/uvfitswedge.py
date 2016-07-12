import time
from astropy.io import fits
import os,sys
import numpy as n
from pylab import *
import capo
import aipy as a
from capo import pspec
import progressbar
poles = ['xx','yy','xy','yx']
WINDOW='blackman-harris'
maxbl = 1000
minbl = 0.01
timeprint=False
windowpad = 0.1 #window padding in Mpc^-1
int_count = 14 # sum this many integrations together before cleaning 
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
    return n.arange(f0-df*fi0,f0-df*fi0 + nfreqs*df,df)
def length(A):
    return n.sqrt(n.dot(A,A))
def gentle(dd, _w, tol=1e-1, area=None, stop_if_div=False,
                    maxiter=100,verbose=False):

    cc, info = a.deconv.clean(dd, _w, tol=tol, area=area, stop_if_div=False,
                    maxiter=100,verbose=False)
    dd = info['res']
    inside_res = n.std(dd[area!=0])
    outside_res = n.std(dd[area==0])
    initial_res = inside_res
    #print inside_res,'->',  
    ncycle=0
    while(inside_res>outside_res):
        _d_cl, info = a.deconv.clean(dd, _w, tol=tol, area=area, stop_if_div=False,
                    maxiter=100,verbose=False)
        res = info['res']
        inside_res = n.std(res[area!=0])
        outside_res = n.std(res[area==0])
        dd = info['res']
        cc += _d_cl
        ncycle += 1
    info['ncycle'] = ncycle
    info['initial_residual'] = initial_res
    info['final_residual'] = inside_res 
    #print inside_res, "N = ",ncycle
    return cc,info
#BEGIN SETTING UP ALL THE STUFF
if len(sys.argv)<2:sys.exit()
F = fits.open(sys.argv[1])
D = F[0]
times = D.data.field('DATE')
t_int = n.diff(times).max()*24*3600
print "t_int = ",t_int
totalt = n.ceil((times.max() - times.min())*24*3600 + t_int)
print "total observing time = ",totalt
Ntimes = n.ceil(totalt/t_int)
print "ntimes = ",Ntimes
bls = D.data.field('BASELINE')
ant2 = (bls%256).astype(int)
ant1 = ((bls-ant2)/256).astype(int)
uvws  = n.array(zip(D.data.field('UU'),D.data.field('VV'),D.data.field('WW')))*1e9 #bl vectors in ns
Nants = len(set(ant2))
Nmax = ant2.max()
freqs = gethdufreqs(D)
z = pspec.f2z(n.mean(freqs)/1e9)
print "padding the horizon by k_parr = ",windowpad
windowpad_delay = windowpad / pspec.dk_deta(z)
print "at z=%4.1f this is %7.5f ns"%(z,windowpad_delay)
df = n.diff(freqs)[0]
delays  =   n.fft.fftfreq(freqs.shape[0],df/1e9) #delays in ns
print "n_ant = ",Nants
print n.sum((ant2-ant1)==0)
Nblt = D.data.field('DATA').shape[0]
Nbl = Nblt/Ntimes
MAXBL = n.max(bls)
Nfreqs = D.data.field('DATA').shape[3]
Npol = D.data.field('DATA').shape[4]
#form up a list of baselines, sorted by length
bl_lengths = []
bls_sorted = []
for i in range(Nblt):
    if bls[i] in bls_sorted:continue
    if length(uvws[i])>maxbl or length(uvws[i])<minbl:continue
    bl_lengths.append(length(uvws[i]))
    bls_sorted.append(bls[i])
#sort everything by baseline length    
I = n.argsort(bl_lengths)
bls_sorted = n.array(bls_sorted)[I]
ant1_sorted = ant1[I]
ant2_sorted = ant2[I]
bl_lengths = n.array(bl_lengths)[I]

#initialize the power spectrum accumulator 
P = n.zeros((len(bl_lengths),Nfreqs))
P_res = n.zeros((len(bl_lengths),Nfreqs))
P2 = n.zeros((len(bl_lengths),Nfreqs)) #this accumulates the square
P2_res = n.zeros((len(bl_lengths),Nfreqs))

C = n.zeros((len(bl_lengths),1)) #this is the integration count
tic = time.time()
initial_rms = []
final_rms = []
ncycle = []
for uvfile in sys.argv[1:]:
    outfile = os.path.basename(uvfile).replace('.uvfits','.p')
    F = fits.open(uvfile)
    D = F[0]
    times = D.data.field('DATE')
    t_int = n.diff(times).max()*24*3600
    print "t_int = ",t_int
    totalt = n.ceil((times.max() - times.min())*24*3600 + t_int)
    print "total observing time = ",totalt
    Ntimes = n.ceil(totalt/t_int)
    print "ntimes = ",Ntimes
    bls = D.data.field('BASELINE')
    ant2 = (bls%256).astype(int)
    ant1 = ((bls-ant2)/256).astype(int)
    uvws  = n.array(zip(D.data.field('UU'),D.data.field('VV'),D.data.field('WW')))*1e9 #bl vectors in ns

    window = None
    waterfall =[]
    accum = n.zeros((Nbl,Nfreqs,Npol)).astype(n.complex64)
    waccum = n.zeros_like(accum)
    naccum = n.zeros((Nbl,Npol))
    bar = progressbar.ProgressBar()
    pol = 0 
    for i in bar(range(Nblt)):
        bl = bls[i]
        bli = n.abs(bls_sorted-bl).argmin() #find the baseline length index
        if ant1[i]==ant2[i]:continue
        #if not (ant1[i]==a2i(30) and ant2[i]==a2i(93)):continue   
        uvw = uvws[i]
        bldelay = length(uvw) + windowpad_delay
	#print "horizon:",length(uvw)," + extra",windowpad_delay
        if bldelay>maxbl or bldelay<minbl:continue #throw out bls longer than 600m
        #grab the data from the big chunk of data in RAM
        data = D.data.field('DATA')[i,:,:,:,pol,0] + 1j*D.data.field('DATA')[i,:,:,:,pol,1] 
        mask = (D.data.field('DATA')[i,:,:,:,pol,2]==0).squeeze()
        data = n.ma.masked_where(mask,data.squeeze())
        if naccum[bli,pol]<int_count: #XXX
            accum[bli,:,pol] += data
            waccum[bli,:,pol] += 1 - data.mask.astype(n.float)
            naccum[bli,pol] += 1
            continue
        else:
            data += accum[bli,:,pol]
            w = waccum[bli,:,pol] + 1 - data.mask.astype(n.float)
            data[w>0] /= w[w>0]
            mask = (naccum[bli,pol]*0.3>w) #mask when we have 30% or less data
            naccum[bli,pol] = 0
        #compute the psf matrix
        #print mask
        _w = n.fft.ifft(1-mask)
        gain = n.abs(_w[0])  
        if gain==0:continue #skip on if theres no data!
        if window is None:
            window = a.dsp.gen_window(data.shape[0],WINDOW)
        dd      =   n.fft.ifft(window*data)
        area = n.zeros_like(dd).astype(n.int)
        area[n.abs(delays)<bldelay] = 1
        #print "Number of in-horizon delay bins",
        #print n.sum(area)
        #_d_cl, info = a.deconv.clean(dd, _w, tol=1e-5, area=area, stop_if_div=False,
        #    maxiter=1000,verbose=False)
        tic = time.time()
        _d_cl, info = gentle(dd, _w, tol=1e-1, area=area, stop_if_div=False,
            maxiter=100,verbose=False)
        toc = time.time()
        if not timeprint: 
            timeprint=True
            print "estimated completion time:",Nblt*(toc-tic)/60,"m"
        initial_rms.append(info['initial_residual'])
        final_rms.append(info['final_residual'])
        ncycle.append(info['ncycle'])
        dd = n.fft.fftshift(dd)
        _d_cl = n.fft.fftshift(_d_cl)
        res = n.fft.fftshift(info['res'])
        try:
            P[bli,:] += _d_cl + res
            P2[bli,:] += (_d_cl+ res) * n.conj(_d_cl+ res)
            P_res[bli,:] += res
            P2_res[bli,:] += res*n.conj(res)
            C[bli] += 1
        except IndexError as e:
            print bl
            raise(e)
        if ant1[i]==a2i(30) and ant2[i]==a2i(93):
            waterfall.append(res)
    #figure()
    waterfall = n.array(waterfall)
    n.savez(outfile+'.30_93.npz',waterfall=waterfall)
    #print waterfall.shape,waterfall.dtype
    #subplot(311)
    #imshow(n.log10(n.abs(waterfall)),aspect='auto')
    #colorbar()
    #subplot(312)
    #plot(waterfall.T)
    #subplot(313)
    #plot(n.mean(waterfall,axis=0))
    #show()
    #save the sum,var and counts for later multiplication
    print outfile+'.npz'
    n.savez(outfile,P=P,P2=P2,P_res=P_res,P2_res=P2_res,C=C,delays=delays,bl_lengths=bl_lengths,freq=n.mean(freqs),ant1=ant1_sorted,ant2=ant2_sorted,initial_rms=n.array(initial_rms),final_rms=n.array(final_rms),ncycle=n.array(ncycle))
#P[C>0] /= C[C>0]
#P2[C>0] /= C[C>0]

#PROD = (P*n.conj(P) - P2)
#PROD[C>0] /= 2*C[C>0]
#PROD = n.sqrt(PROD)
#
#
print "time to scan bls and dt = ",time.time()-tic
#imshow(n.log10(n.abs(PROD.T)),aspect='auto')
#n.savez('mywedge',P=PROD,delays=delays,bl_lengths=bl_lengths,freq=n.mean(freqs))
#colorbar()
#show()
