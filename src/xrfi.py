'''Module for all things Radio Frequency Interference Flagging'''
import numpy as np

def medmin(chisq):
    #return np.median(np.min(chisq,axis=0))
    mn = np.min(chisq,axis=0)
    return 2*np.median(mn) - np.min(mn)

def omni_chisq_to_flags(chisq, K=8, sigma=6, sigl=2):
    '''Returns a mask of RFI given omnical's chisq statistic'''
    w_sm = np.empty_like(chisq)
    sig = np.empty_like(chisq)
    #get smooth component of chisq
    for i in xrange(chisq.shape[0]):
        for j in xrange(chisq.shape[1]):
            i0,j0 = max(0,i-K), max(0,j-K)
            i1,j1 = min(chisq.shape[0], i+K), min(chisq.shape[1], j+K)
            #w_sm[i,j] = np.median(chisq[i0:i1,j0:j1])
            w_sm[i,j] = medmin(chisq[i0:i1,j0:j1])
    #the residual from smooth component
    w_rs = chisq - w_sm 
    w_sq = np.abs(w_rs)**2
    #get the standard deviation of the media.
    for i in xrange(chisq.shape[0]):
        for j in xrange(chisq.shape[1]):
            i0,j0 = max(0,i-K), max(0,j-K)
            i1,j1 = min(chisq.shape[0], i+K), min(chisq.shape[1], j+K)
            #sig[i,j] = np.sqrt(np.median(w_sq[i0:i1,j0:j1]))
            sig[i,j] = np.sqrt(medmin(w_sq[i0:i1,j0:j1]))
    #Number of sigma above the residual unsmooth part is.
    f1 = w_rs / sig
    #mask off any points above 'sig' sigma and nan's.
    f1 = np.ma.array(f1, mask=np.where(f1>sigma,1,0)) 
    f1.mask | np.isnan(f1)
    
    #Start the watershed
    prevx = 0
    prevy = 0
    x,y = np.where(f1.mask)
    while x.size != prevx and y.size != prevy:
        for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
            prevx = x.size
            prevy = y.size
            xp, yp = (x+dx).clip(0,f1.shape[0]-1), (y+dy).clip(0,f1.shape[1]-1)
            i = np.where(f1[xp,yp]>sigl)[0] # if sigma > 'sigl'
            f1.mask[xp[i],yp[i]] = 1
            x,y = np.where(f1.mask)
    
    f1ch = np.average(f1.mask, axis=0); f1ch.shape = (1,-1)
    #The cut off value is a made up number here...sig = 'sig' if none flagged.
    f1.mask = np.logical_or(f1.mask, np.where(f1 > sigma*(1-f1ch), 1, 0))
    f1t = np.average(f1.mask, axis=1)
    ts = np.where(f1t > 2*np.median(f1t))
    f1.mask[ts] = 1
    f1f_sum = np.sum(f1.filled(0), axis=0)
    f1f_wgt = np.sum(np.logical_not(f1.mask), axis=0)
    f1f = f1f_sum / f1f_wgt.clip(1,np.Inf)
    fs = np.where(f1f > 2)
    f1.mask[:,fs] = 1
    mask = f1.mask
    return mask
