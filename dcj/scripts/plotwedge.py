from matplotlib import use
use('Agg')
from pylab import *
import numpy as n
import time
from capo import pspec
from astropy.time import Time


#STUFF RELATED TO PUTTING INFO ON THE PLOTS
import mwapy
from mwapy import ephem_utils
from mwapy.get_observation_info import *
from mwapy.obssched.base import schedule

# configure the logging
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('get_observation_info')
logger.setLevel(logging.WARNING)

# open up database connection
try:
    db = schedule.getdb()
except:
    logger.error("Unable to open connection to database")
    sys.exit(1)


rcParams['font.size'] = 22
rcParams['font.family']='serif'
NBLS  = 2000
tic = time.time()
F = n.load(sys.argv[1])
bl_lengths = F['bl_lengths']
Nbl  = len(bl_lengths)
delays = F['delays']
freq = F['freq']
windowpad = F['windowpad']
plotres=True
PLOTWINDOWCONTOUR=False
figsize=(20,25)
try: 
    observation_num  = int(sys.argv[1].split('.')[0])
    observation=MWA_Observation(observation_num,db=db)

except(ValueError):
    gpstime = int(n.round(Time(F['time'],scale='utc',format='jd').gps))
    observation_num = find_closest_observation(gpstime,db=db)
    observation=MWA_Observation(observation_num,db=db)


for plotres in [False,True]:
    z = pspec.f2z(freq/1e9)
    #windowpad = 0.01 #should match whatever I used in uvfitswedge.py
    if len(sys.argv[1:])==1 and not sys.argv[1].endswith('.p.npz'):
        amps = n.ma.masked_invalid(n.log10(n.abs(F['P'].T)))
        outname = sys.argv[1].replace('.npz','')
    else:
        if len(sys.argv[1:])==1:outname = sys.argv[1].replace('.npz','')
        else:outname = 'mywedge'
        if plotres: outname += '_res'
        P=0
        P2 = 0
        C = 0
        for Pfile in sys.argv[1:]:
            F = n.load(Pfile)
    	    if plotres:
    	        P += F['P_res']
                P2 += F['P2_res']
            else:
                P += F['P']
                P2 += F['P2']
            C += F['C']        
        PROD = P*n.conj(P) - P2
        #correctshape = PROD.shape
        #PROD.shape = C.shape+(PROD.shape[1],)
        #print PROD[C>0].shape
        PROD[C>0] /= 2*C[C>0]


        #PROD.shape = correctshape
        PROD = n.swapaxes(n.sqrt(n.abs(PROD)),1,0)

        P = PROD.copy()
        amps = n.ma.masked_invalid(n.log10(PROD))
        windowpad_ns = F['windowpad']
    B,D = n.meshgrid(bl_lengths,n.fft.fftshift(delays))
    print B.shape,amps.shape
    print "nancount", n.sum(n.isnan(P)),'post div'

    if True:
        fig = figure(figsize=figsize)
        ax1 = subplot(211)
        title('XX')
        pcolor(B,D,amps[:,:,0],vmin=0)
        plot(bl_lengths,bl_lengths+windowpad,'k')
        plot(bl_lengths,-1*bl_lengths-windowpad,'k')
        ylim([-2*bl_lengths.max(),2*bl_lengths.max()])
        ylabel('delay [ns]')
        xlabel('bl length [ns]')
        colorbar()
        ax2 = subplot(212)
        title('YY')
        pcolor(B,D,amps[:,:,1],vmin=0)
        plot(bl_lengths,bl_lengths+windowpad,'k')
        plot(bl_lengths,-1*bl_lengths-windowpad,'k')
        ylim([-2*bl_lengths.max(),2*bl_lengths.max()])
        ylabel('delay [ns]')
        xlabel('bl length [ns]')
        colorbar()
        print "time = ",time.time()-tic
        savefig(outname+'.dwedge.png')
        suptitle(sys.argv[1].split('.')[0] + \
        ', LST=' + str(n.round(observation.LST,2))+', HA=' + str(n.round(observation.HA,2))+ \
        ', AZ/EL='+str(n.round(observation.azimuth,2))+\
        '/'+str(n.round(observation.elevation,2)))
        savefig(outname+'.dwedge_info.png')
    
    #average +/- k//
    #plot with k units
    #add bandwidth to horizon limit
    figure(figsize=figsize)
    ndelays = len(delays)
    nkparr = n.sum(delays>=0)
    #P = F['P'].T
    #FOLD UP THE PSPEC AROUND ZERO RELATIVE DELAY!!!
    #XXX THIS WILL BE INCORRECT FOR PHASED UP STUFF.
    wedge_X = n.sqrt(n.abs(P[nkparr:,:,0]*n.conj(P[nkparr:,:,0]) + \
        n.flipud(P[:nkparr,:,0]*n.conj(P[:nkparr,:,0]))/2))
    wedge_Y = n.sqrt(n.abs(P[nkparr:,:,1]*n.conj(P[nkparr:,:,1]) + \
        n.flipud(P[:nkparr,:,1]*n.conj(P[:nkparr,:,1]))/2))
    print n.sum(n.isnan(P)),n.sum(P<0)
    
    wedge_X = n.ma.masked_invalid(n.log10(wedge_X))
    wedge_Y = n.ma.masked_invalid(n.log10(wedge_Y))
    
    
    kparr = delays[delays>=0]*pspec.dk_deta(z)
    kperp = bl_lengths*pspec.dk_du(z)
    horizon = bl_lengths*pspec.dk_deta(z)
    print bl_lengths,pspec.dk_deta(z),z
    print horizon
    windowpad = windowpad_ns*pspec.dk_deta(z)
    KPERP,KPARR = n.meshgrid(kperp,kparr)
    MIN = n.min([n.mean(wedge_X[KPARR>(horizon+windowpad)]),n.mean(wedge_X[KPARR>(horizon+windowpad)])])*0.8
    
    dfcoarse=30/24.*1e-3
    if PLOTWINDOWCONTOUR:
        #calculate danny's stat window
        BL,D = n.meshgrid(bl_lengths,delays[delays>=0])
        coarse_chan_delay =  1/(1.28/1e3)#
        print "coarse_chan_delay",coarse_chan_delay
        WINDOW = n.zeros_like(D)
        WINDOW[n.logical_and(D>(windowpad_ns + bl_lengths),D<coarse_chan_delay)] = 1
        print "window bins included:",n.sum(WINDOW)


    fig = figure(figsize=figsize)
    ax1 = subplot(211)
    title('XX')
    pcolor(KPERP,KPARR,wedge_X,vmin=MIN)
    if PLOTWINDOWCONTOUR:
        contour(KPERP,KPARR,WINDOW,[0.9],linewidths=[4],colors='r')
    ylim([0,kparr.max()/2])
    xlim([kperp.min(),kperp.max()])
    plot(kperp,horizon,'k')
    ylabel('$k_\\parallel$ [hMpc$^{-1}$]')
    xlabel('$k_\\perp$ [hMpc$^{-1}$]')
    twiny()
    xlim([bl_lengths.min(),bl_lengths.max()])
    xlabel('baseline length [$\lambda$]')
    ax2 = subplot(212)
    title('YY')
    pcolor(KPERP,KPARR,wedge_Y,vmin=MIN)
    if PLOTWINDOWCONTOUR:
        contour(KPERP,KPARR,WINDOW,[0.9],linewidths=[4],colors='r')
    ylim([0,kparr.max()/2])
    xlim([kperp.min(),kperp.max()])
    plot(kperp,horizon,'k')
    ylabel('$k_\\parallel$ [hMpc$^{-1}$]')
    xlabel('$k_\\perp$ [hMpc$^{-1}$]')
    twiny()
    xlim([bl_lengths.min(),bl_lengths.max()])
    xlabel('baseline length [$\lambda$]')
    
    
    
    savefig(outname+'.wedge.png')
        


    print "time = ",time.time()-tic
    #show()
