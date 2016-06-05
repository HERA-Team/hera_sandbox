from matplotlib import use
use('Agg')
from pylab import *
import numpy as n
import time,sys
from capo import pspec
rcParams['font.size'] = 22
rcParams['font.family']='serif'
SHORTWEDGEBL_LENGTH  = 50 #lambda
tic = time.time()


def fold_wedge(P,nkparr):
    #fold P in the first dimension at index nkparr
    return n.sqrt(n.abs(P[nkparr:,:] + \
    n.flipud(n.conj(P[:nkparr,:]))/2))
def p(x):
    return n.round(x,5)
#load the wedge data
for Pfile in sys.argv[1:]:

    P=0
    P2 = 0
    P_res = 0
    P2_res = 0
    C = 0
    F = n.load(Pfile)
    bl_lengths = F['bl_lengths']
    delays = F['delays']
    freq = F['freq']
    z = pspec.f2z(freq/1e9)
    NBLS=len(F['bl_lengths'])
    windowpad_ns = F['windowpad']
    P_res += F['P_res']
    P2_res += F['P2_res']
    P += F['P']
    P2 += F['P2']
    C += F['C']        
    F.close()
    PROD = P*n.conj(P) - P2
    PROD_res = P_res*n.conj(P_res) - P2_res
    PROD[C>0] /= 2*C[C>0]
    PROD_res[C>0] /= 2*C[C>0]
    PROD = n.swapaxes(PROD,1,0)
    PROD_res = n.swapaxes(PROD_res,1,0)
    
    ndelays = len(delays)
    nkparr = n.sum(delays>=0)
    #FOLD UP THE PSPEC AROUND ZERO RELATIVE DELAY!!!
    #XXX THIS WILL BE INCORRECT FOR PHASED UP STUFF.
    wedge_X = fold_wedge(PROD[:,:,0],nkparr)
    wedge_X_res = fold_wedge(PROD_res[:,:,0],nkparr)
    wedge_Y= fold_wedge(PROD[:,:,1],nkparr)
    wedge_Y_res= fold_wedge(PROD_res[:,:,1],nkparr)
    #wedge_X = n.sqrt(n.abs(P[nkparr:,:,0]*n.conj(P[nkparr:,:,0]) + \
    #    n.flipud(P[:nkparr,:,0]*n.conj(P[:nkparr,:,0]))/2))
    #wedge_Y = n.sqrt(n.abs(P[nkparr:,:,1]*n.conj(P[nkparr:,:,1]) + \
    #    n.flipud(P[:nkparr,:,1]*n.conj(P[:nkparr,:,1]))/2))
    #print n.sum(n.isnan(P)),n.sum(P<0)
    
    #wedge_X = n.ma.masked_invalid(n.log10(wedge_X))
    #wedge_Y = n.ma.masked_invalid(n.log10(wedge_Y))
    
    kparr = delays[delays>=0]*pspec.dk_deta(z)
    BL,D = n.meshgrid(bl_lengths,delays[delays>=0])
    kperp = bl_lengths*pspec.dk_du(z)
    horizon = bl_lengths*pspec.dk_deta(z)
    #print bl_lengths,pspec.dk_deta(z),z
    #print horizon
    windowpad = windowpad_ns*pspec.dk_deta(z)
    coarse_chan_delay =  1/(1.28/1e3)#
    WINDOW = n.zeros_like(D)
    WINDOW[n.logical_and(D>(windowpad_ns + bl_lengths),D<coarse_chan_delay)] = 1
    n.set_printoptions(precision=2)
    WEDGE = n.zeros_like(D)
    WEDGE[D<(windowpad_ns + bl_lengths)] = 1
    SHORTWEDGE = n.zeros_like(D)
    SHORTWEDGE[n.logical_and(D<(windowpad_ns + bl_lengths),BL<50)] = 1
    LONGWEDGE = n.zeros_like(D)
    LONGWEDGE[n.logical_and(D<(windowpad_ns + bl_lengths),BL>50)] = 1
    

    print Pfile.split('.')[0], #obsid
    print p(n.mean(wedge_X[WEDGE>0])),p(n.mean(wedge_Y[WEDGE>0])), #WEDGE POWR
    print p(n.mean(wedge_X_res[WINDOW>0])),p(n.mean(wedge_Y_res[WINDOW>0])), #WINDOW POWER
    print p(n.mean(wedge_X_res[WEDGE>0])),p(n.mean(wedge_Y_res[WEDGE>0])),  #RESIDUAL WEDGE 
    print p(n.mean(wedge_X_res[SHORTWEDGE>0])),p(n.mean(wedge_Y_res[SHORTWEDGE>0])), #GALACTIC WEDGE
    print p(n.mean(wedge_X_res[LONGWEDGE>0])),p(n.mean(wedge_Y_res[LONGWEDGE>0]))    #POINT SOURCE WEDGE
    sys.stdout.flush()

#KPERP,KPARR = n.meshgrid(kperp,kparr)
#MIN = n.min([n.mean(wedge_X[KPARR>(horizon+windowpad)]),n.mean(wedge_X[KPARR>(horizon+windowpad)])])*0.8
#average the power above the clean box (horizon + padding) and below the first 1.28MHz harmonic

