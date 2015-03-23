from matplotlib import use
use('Agg')
import time
from astropy.io import fits
import os,sys
import numpy as n
import capo
import aipy as a
from capo import pspec
import progressbar
from scipy.io import readsav
import ephem
import resource
from glob import glob
poles = ['xx','yy','xy','yx']
WINDOW='blackman-harris'
maxbl = 1000
minbl = 0.01
ns_per_m = 3.33564095
singlebl = False
VERB=False
VERBCLEAN=False
timeprint=False
windowpad = 0.01 #window padding in Mpc^-1
windowpad_ns = 60  # any additional window padding in ns
int_count = 14 # sum this many integrations together before cleaning 
dirty = False
dophase=True
T1 = 18
T2 = 48
def ram():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024 #return total ram usage in GB
class MWA_pos():
    def __init__(self,antposfile='MWA_128T_antenna_locations.txt'):
        self.antpos = n.loadtxt(antposfile)*ns_per_m
    def get_index_uvw(self,i,j):
        #input mwa ant numbers in 0 index
        #returns baseline uvw in ns 
        return self.antpos[i,:] - self.antpos[j,:]
    def get_antnum_uvw(self,ant1,ant2):
        #input mwa ant numbers 11-168
        #returns baseline uvw vector in ns
        i,j = a2i(ant1),a2i(ant2)
        return self.get_index_uvw(i,j)

if singlebl:
    from pylab import *
def bodyxyz(body):
    x = n.cos(body.alt)*n.cos(body.az)
    y = n.cos(body.alt)*n.sin(body.az)
    z = n.sin(body.alt)
    return n.array([x,y,z])
def calcdelays(phase_center,uvws,times):
    print "phasing to xyz ",
    delays = []
    O = ephem.Observer()
    O.lat = '-26:42:11.95'
    O.lon = '116:40:14.93'
    S = ephem.FixedBody()
    S._ra = phase_center[0] * n.pi/180
    S._dec = phase_center[1]*n.pi/180
    O.date = times[0] - 2415020 #convert from julian to dublin julian
    S.compute(O)
    #print S._dec
    #print S.dec
    #print S._ra
    #print S.ra
    ##print S.ra,S.dec
    #print S.alt,S.az
    print bodyxyz(S)
    for i,uvw in enumerate(uvws):
        O.date = times[0] - 2415020
        S.compute(O)
        s = bodyxyz(S)
        delays.append(n.dot(uvw,s))        
    delays = n.array(delays)
    return delays    
def i2a(i):
    #convert from MWA to unit ant index
    # assumes zero indexed antenna numbers 0-127
    tens = int(i/8)+1
    ones = i%8+1
    return tens*10+ones
def a2i(a):
    #convert from unit to MWA ant index
    #returns a zero indexed number 0-127
    eights = int(a/10)-1
    ones = a%10
    return eights*8+ones-1
i2a = n.vectorize(i2a)
a2i = n.vectorize(a2i)
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
    if verbose:
        print "entering gentle clean"
    cc, info = a.deconv.clean(dd, _w, tol=tol, area=area, stop_if_div=False,
                    maxiter=maxiter,verbose=verbose)
    cc = n.zeros_like(dd)
    inside_res = n.std(dd[area!=0])
    outside_res = n.std(dd[area==0])
    initial_res = inside_res
    #print inside_res,'->',  
    ncycle=0
    if verbose:
        print "inside_res outside_res"
        print inside_res, outside_res
    inside_res = 2*outside_res #just artifically bump up the inside res so the loop runs at least once
    while(inside_res>outside_res and maxiter>0):
        if verbose: print '.',
        _d_cl, info = a.deconv.clean(dd, _w, tol=tol, area=area, stop_if_div=False,
                    maxiter=maxiter,verbose=verbose)
        res = info['res']
        inside_res = n.std(res[area!=0])
        outside_res = n.std(res[area==0])
        dd = info['res']
        cc += _d_cl
        ncycle += 1
        if verbose: print inside_res,outside_res
        if ncycle>1000:break
    info['ncycle'] = ncycle
    info['initial_residual'] = initial_res
    info['final_residual'] = inside_res 
    return cc,info
def loadfhdsav(savfile):
    #get the other poles
    pol = savfile.split('_')[-1][:2]
    if pol=='xx':
        uvfile_altpol = savfile.replace('xx','yy')
        pols = [0,1]
    elif pol=='yy':
        uvfile_altpol = savfile.replace('yy','xx')
        pols = [1,0]
    else:
        print "polarization not found in filename. skipping"
        raise(StandardError)
    if not os.path.exists(uvfile_altpol):
        print "pol file",uvfile_altpol,"not found. please find"
        raise(StandardError)
    #paramfile = savfile.split('_')[0]+'_params.sav'
    paramfile = savfile.replace('vis_%s'%(['xx','yy'][pols[0]]),'params')
    if not os.path.exists(paramfile):
        print "error: paramfile=",paramfile,"not found. please find"
        raise(StandardError)
    #weightfile = savfile.split('_')[0]+'_flags.sav'
    weightfile = savfile.replace('vis_%s'%(['xx','yy'][pols[0]]),'flags')
    if not os.path.exists(weightfile):
        print "error: weightfile",weightfile,"not found, please find"
        raise(StandardError)
    print "loading:",savfile
    uvfile = readsav(savfile)
    ant1 =  uvfile['obs']['baseline_info'][0]['tile_a'][0]-1
    print "min(ant1)=",ant1.min(),"should be 0"
    ant2 =  uvfile['obs']['baseline_info'][0]['tile_b'][0]-1
    print "max(ant2)=",ant2.max(),"should be 127"
    data =  uvfile['vis_ptr']
    #times = uvfile['obs']['baseline_info'][0]['jdate'][0]
    baselines = (ant2)*256+ant1
    freqs = uvfile['obs']['baseline_info'][0]['freq'][0]

    print "loading alternate polarization",uvfile_altpol
    uv_altpol = readsav(uvfile_altpol)
    data_altpol = uv_altpol['vis_ptr']

    print "loading baselines from params file:",paramfile
    params = readsav(paramfile)
    U = params['params']['uu'][0]*1e9
    V = params['params']['vv'][0]*1e9
    W = params['params']['ww'][0]*1e9
    uvw = n.array(zip(U,V,W))
    times = params['params']['time'][0]
    print "loading weights from :",weightfile
    flags = readsav(weightfile)
    mask = n.dstack([flags['flag_arr'][0],flags['flag_arr'][1]])==0 #the zeros are the flags 
    #create the new fits file
    outdata = n.zeros((data.shape[0],data.shape[1],2)).astype(n.complex64)
    outdata[:,:,pols[0]] = data
    outdata[:,:,pols[1]] = data_altpol
    ant1,ant2 = i2a(ant1),i2a(ant2)
    return uvw,ant1,ant2,baselines,times,freqs,outdata,mask
#BEGIN SETTING UP ALL THE STUFF 
#  COMMENTED BECAUSE I MOVED ALL THIS INTO THE FILE LOOP.
#  ToBeDELETED
#    LOAD THE DATA INTO RAM
#if len(sys.argv)<2:sys.exit()
#if sys.argv[1].endswith('.uvfits'):
#    print "initializing variables with ",sys.argv[1]
#    F = fits.open(sys.argv[1])
#    D = F[0]
#    Nblt = D.data.field('DATA').shape[0] 
#    times = D.data.field('DATE')
#    uvws  = n.array(zip(D.data.field('UU'),D.data.field('VV'),D.data.field('WW')))*1e9 #bl vectors in ns
#    bls = D.data.field('BASELINE')
#    freqs = gethdufreqs(D)
#    Nfreqs = D.data.field('DATA').shape[3]
#    Npol = D.data.field('DATA').shape[4]
#    phase_center = [D.header['CRVAL5'],D.header['CRVAL6']]
#elif sys.argv[1].endswith('.sav'):
#    uvws,ant1,ant2,bls,times,freqs,data,weights = loadfhdsav(sys.argv[1])
#    Npol = 2
#    Nfreqs = len(freqs)
#    Nblt = data.shape[0]
#    phase_center = None
if not dophase:
    phase_center=None
print "total ram usage = ",ram(),"GB"

MWA_baselines = MWA_pos('/nfs/blank/h4215/djacobs/MWA_128T_antenna_locations.txt')#load the antenna locations for computing baselines

for uvfile in sys.argv[1:]:
    print uvfile, 
    if not os.path.exists(uvfile):
        print "file not found!"
        continue
    print "-->",
    #build the name of the output file
    if uvfile.endswith('uvfits'):
        outfile = os.path.basename(uvfile).replace('.uvfits','.p')
    else:
        outfile = os.path.basename(uvfile).split('_')[0]+'.fhd.p'
    if dirty: outfile +='.dirty'
    print outfile+'.npz'
    sys.stdout.flush()
    if os.path.exists(outfile):
        print "... exists. Skipping."
        continue
    #LOAD THE DATA
    if uvfile.endswith('uvfits'):
        F = fits.open(uvfile)
        D = F[0]
        times = D.data.field('DATE')
        bls = D.data.field('BASELINE')
        uvws  = n.array(zip(D.data.field('UU'),D.data.field('VV'),D.data.field('WW')))*1e9 #bl vectors in ns
        DATA =D.data.field('DATA').squeeze()[:,:,:,0] + 1j*D.data.field('DATA').squeeze()[:,:,:,1]
        MASK = (D.data.field('DATA').squeeze()[:,:,:,2]==0)
        freqs = gethdufreqs(D)
        Nfreqs = D.data.field('DATA').shape[3]
        Npol = D.data.field('DATA').shape[4]
        Nblt = D.data.field('DATA').shape[0] 
        del(D)
        ant2    = n.array(bls%256).astype(int)-1
        ant1    = n.array((bls-ant2)/256).astype(int)-1
        ant1,ant2 = n.array(map(i2a,ant1)),n.array(map(i2a,ant2))
    else:
        uvws,ant1,ant2,bls,times,freqs,DATA,MASK = loadfhdsav(uvfile)
        print "fraction of zeros in MASK",n.sum(MASK==0)/float(MASK.size)
    print "useful things about antenna numbering"
    print "the largest antenna number is ",ant2.max()
    print "the smalleest antenna number is",ant1.min()
    print "the range of the MASK variable is", MASK.min(),MASK.max()
    Nblt = DATA.shape[0]
    print "Nblt = ",Nblt,len(times)
    Npol = DATA.shape[-1]
    Nfreqs = DATA.shape[-2]
    print "finished loading data, deleting intermediate object"
    
    print "total ram usage = ",ram(),"GB"

    """
    Setup a bunch of variables
    """
    #### BEGIN VARIABLE SETUP
    t_int   = n.diff(times).max()*24*3600
    print "t_int = ",t_int
    totalt  = n.ceil((times.max() - times.min())*24*3600 + t_int)
    print "total observing time = ",totalt
    Ntimes  = n.ceil(totalt/t_int)
    Nants   = len(set(ant2))
    Nmax    = ant2.max()
    z       = pspec.f2z(n.mean(freqs)/1e9)
    print "padding the horizon by k_parr = ",windowpad+windowpad_ns
    windowpad_delay = windowpad / pspec.dk_deta(z) + windowpad_ns
    print "at z=%4.1f this is %7.5f ns"%(z,windowpad_delay)
    df      = n.diff(freqs)[0]
    delays  =   n.fft.fftfreq(freqs.shape[0],df/1e9) #delays in ns
    print "n_ant = ",Nants
    print "number of auto correlations = ", n.sum((ant2-ant1)==0)
    Nbl = Nblt/Ntimes
    MAXBL = n.max(bls)

    #form up a list of baselines, sorted by length
    bl_lengths = []
    bls_sorted = []
    uvw_sorted = []
    ant1_sorted = []
    ant2_sorted = []
    blcount = {}
    for i in range(Nblt):
        if bls[i] in bls_sorted:continue
        if length(uvws[i])>maxbl or length(uvws[i])<minbl:continue
        #if ant1[i]==28 or ant2[i]==28:
        #    print ant1[i],ant2[i],uvws[i],bls[i],minbl,maxbl    
        bl_lengths.append(length(uvws[i]))
        bls_sorted.append(bls[i])
        uvw_sorted.append(uvws[i])
        ant1_sorted.append(ant1[i])
        ant2_sorted.append(ant2[i])
    #sort everything by baseline length    
    I = n.argsort(bl_lengths)
    bls_sorted = n.array(bls_sorted)[I]
    bl_lengths = n.array(bl_lengths)[I]
    uvw_sorted = n.array(uvw_sorted)[I]
    ant1_sorted = n.array(ant1_sorted)[I]
    ant2_sorted = n.array(ant2_sorted)[I]
    ##### END VARIABLE SETUP
  
    #### BEGIN POWER SPECTRUM COMPUTATION
    #initialize the power spectrum accumulator 
    P       =   n.zeros((len(bl_lengths),Nfreqs,Npol)).astype(n.complex64) 
    P_res   =   n.zeros((len(bl_lengths),Nfreqs,Npol)).astype(n.complex64) 
    P2      =   n.zeros((len(bl_lengths),Nfreqs,Npol)).astype(n.complex64)  #this accumulates the square
    P2_res  =   n.zeros((len(bl_lengths),Nfreqs,Npol)).astype(n.complex64) 
    C = n.zeros((len(bl_lengths),1,Npol)) #this is the integration count
    tic = time.time()
    initial_rms = []
    final_rms = []
    ncycle = []




    #CALC the phase delays
#    print "calculating the delays"
#    if not phase_center is None:
#        phase_delays = calcdelays(phase_center,uvws,times)
#    else:
#        phase_delays = n.zeros_like(times)
    if singlebl and VERB: print "plotting a single bl. found nrecords= ",n.sum(n.logical_and(ant1==T1,ant2==T2))
    #CALC THE CONJS
    conj = []
    print "calculating the conjugations"
    for blt,(a1,a2)  in enumerate(zip(ant1,ant2)):
        uvw = MWA_baselines.get_antnum_uvw(a2,a1)
        if VERB and a1==T1 and a2==T2: print "file uvw = ",uvws[blt],"calculated uvw = ",uvw
        ang = n.angle(uvw[0]+1j*uvw[1],deg=True)
        if (ang<(90+22.5)) & (ang>(-90+22.5)):
            conj.append(0)
            #uvws[blt] = uvw
        else:
            conj.append(1)
            #uvws[blt] = -1*uvw
            if VERB and a1==T1 and a2==T2: print "conjugating with angle",ang
#    for uvw in uvws:
#        #if uvw[0]>0:conj.append(1)
#        ang = n.angle(uvw[0]+1j*uvw[1])*180/n.pi
#        if (ang<(90+22.5)) & (ang>(-90+22.5)):conj.append(1)
#        else:conj.append(0)
    conj = n.array(conj)



    
    #INITIALIZE THE ACCUMULATORS!
    window = None
    waterfall_res =[]
    waterfall_cc = []
    accum = n.zeros((Nbl,Nfreqs,Npol)).astype(n.complex64)
    waccum = n.zeros_like(accum)
    naccum = n.zeros((Nbl,Npol))
    #ant1_sorted = n.zeros_like(bls_sorted)
    #ant2_sorted = n.zeros_like(bls_sorted)

    #SCAN THROUGH BASELINES
    bar = progressbar.ProgressBar()
    for i in bar(range(Nblt)):
        bl = bls[i]
        bli = n.abs(bls_sorted-bl).argmin() #find the baseline length index
        if ant1[i]==ant2[i]:continue
        if not (ant1[i]==T1 and ant2[i]==T2) and singlebl:continue   
        if singlebl and VERB: print '.'
        uvw = uvws[i]
        #if singlebl: print "|u| = ",length(uvw)
        bldelay = length(uvw) + windowpad_delay
	#print "horizon:",length(uvw)," + extra",windowpad_delay
        if bldelay>maxbl or bldelay<minbl:
            if singlebl: print "baseline length of ",bldelay," [ns] outside of limit",minbl,"-",maxbl
            continue #throw out bls longer than 600m
        if singlebl and VERB: print '-'
        #SCAN BY POL
        for pol in range(Npol):
            #grab the data from the big chunk of data in RAM
            data = DATA[i,:,pol].squeeze()#D.data.field('DATA')[i,:,:,:,pol,0] + 1j*D.data.field('DATA')[i,:,:,:,pol,1
            if dophase:
                #print freqs
                delay = uvw[2] #the w term points at the current phase center 
                #print delay
                
                data *= n.exp(2*n.pi * 1j * freqs/1e9 * delay)
                if conj[i]: 
                    data = n.conj(data)
                    uvws[bli] *= -1
                    uvw *= -1
            mask = MASK[i,:,pol]#(D.data.field('DATA')[i,:,:,:,pol,2]==0).squeeze()
            if naccum[bli,pol]<(int_count-1): #if we are within the integration window, accumulate
                accum[bli,:,pol] += data
                waccum[bli,:,pol] += 1 - mask.astype(n.float)
                naccum[bli,pol] += 1
                continue
            if naccum[bli,pol]==(int_count-1):
                data += accum[bli,:,pol]
                w = waccum[bli,:,pol] + 1 - mask.astype(n.float)
                data[w>0] /= w[w>0]
                mask = (naccum[bli,pol]*0.3)>w.clip(0,1e10) #mask when we have 30% or less data
                data = n.ma.masked_where(mask,data)
                naccum[bli,pol]    = 0
                accum[bli,:,pol]  *= 0
                waccum[bli,:,pol] *= 0 
                
            else:continue
            if VERB: print "FFT", ant1[bli],ant2[bli], bl,pol,
            #compute the psf matrix
            _w = n.fft.ifft((1-mask.astype(float)))
            if VERB and singlebl: print "nflags = ",n.sum(1-w)
            if VERB and singlebl: print "peak of weight function = ",_w.max()
            gain = n.abs(_w[0])  
            if gain==0 or n.ma.sum(data)==0:
                if VERB:print "skipping due to lack of data",gain, n.ma.sum(data)
                continue #skip on if theres no data!
            elif VERB: print 
            if window is None:
                window = a.dsp.gen_window(data.shape[0],WINDOW)
            dd      =   n.fft.ifft(window*data)  #THE DELAY TRANSFORM
            area = n.zeros_like(dd).astype(n.int)
            area[n.abs(delays)<bldelay] = 1  #THE CLEAN BOX
            #BEGIN CLEANING
            if singlebl: print "CLEANING"
            tic = time.time()
            maxiter=100
            if dirty: maxiter=0
            try:
                _d_cl, info = gentle(dd, _w, tol=1e-1, area=area, stop_if_div=False,
                    maxiter=maxiter,verbose=((singlebl or VERB) and VERBCLEAN))
            except(KeyboardInterrupt):
                print bl,ant1[i],ant2[i]
                _d_cl, info = gentle(dd, _w, tol=1e-1, area=area, stop_if_div=False,
                    maxiter=100,verbose=True)
                sys.exit()
                #PROVIDE A SAFE WAY TO BUST OUT OF STUCK CLEANS
            toc = time.time()
            #END CLEANING. ooooh yeah, we are done!
            if not timeprint: 
                timeprint=True
                print "estimated completion time:",Nblt*(toc-tic)/60,"m"
            initial_rms.append(info['initial_residual'])
            final_rms.append(info['final_residual'])
            ncycle.append(info['ncycle'])
            #if info['ncycle']==maxiter: print "!ncycle!"
            if singlebl and VERB:
                print "finished gentle cleaning"
                print "ncycle = ",info['ncycle']
                print "clean bins = ",n.sum(area)
                print "clean range = ",delays[area>0].min(),delays[area>0].max()
                print "initial residual",initial_rms[-1]
                print "final residual", final_rms[-1]
            dd = n.fft.fftshift(dd)
            _d_cl = n.fft.fftshift(_d_cl)
            res = n.fft.fftshift(info['res'])
            psf = n.fft.fftshift(_w)
            try:
                P[bli,:,pol] += _d_cl + res
                P2[bli,:,pol] += (_d_cl+ res) * n.conj(_d_cl+ res)
                P_res[bli,:,pol] += res
                P2_res[bli,:,pol] += res*n.conj(res)
                C[bli,0,pol] += 1
            except IndexError as e:
                print bli,pol,P.shape
                raise(e)
            #SAVE A WATERFALL FOR PLOTTING PURPOSES
            if ant1[i]==T1 and ant2[i]==T2:
                waterfall_res.append(res)
                waterfall_cc.append(_d_cl)
                if VERB: print 'appending to waterfall npoints:',len(res)
    if len(waterfall_res)>0:
        print "saving waterfall"
        print "delay[0] = ",delays[0]
        waterfall_res = n.array(waterfall_res)
        waterfall_res.shape = (waterfall_res.shape[0]/Npol,Npol,waterfall_res.shape[1])
        waterfall_cc = n.array(waterfall_cc)
        waterfall_cc.shape = (waterfall_cc.shape[0]/Npol,Npol,waterfall_cc.shape[1])

        
        n.savez(outfile+'.mybl.%i_%i.npz'%(T1,T2),waterfall_res=waterfall_res,delays=n.fft.fftshift(delays),waterfall_cc=waterfall_cc,psf=psf)
        print "saving delays too!"
        print "saved ",outfile+'.mybl.%i_%i.npz'%(T1,T2)  
    if VERB: print "total baselines skipped:",n.sum(C==0)
    if not singlebl:
        n.savez(outfile,P=P,P2=P2,P_res=P_res,P2_res=P2_res,C=C,delays=delays,bl_lengths=bl_lengths,freq=n.mean(freqs),ant1=ant1_sorted,ant2=ant2_sorted,initial_rms=n.array(initial_rms),final_rms=n.array(final_rms),ncycle=n.array(ncycle),uvws=uvw_sorted,windowpad=windowpad_delay)
#
    print "time to scan bls and dt = ",time.time()-tic
