#! /usr/bin/env python
'''
Implements a Fringe Rate Filter useing and FIR filter.
'''
import aipy as a
import capo
import numpy as n, pylab as p
import scipy.interpolate
import optparse,sys

def sky_fng_thresh(bl_ew_len, inttime, nints, freq, min_fr_frac=-.3, xtalk=-1,
                        max_fr_frac=1.):
    '''
    For bl_ew_len (the east/west projection) in ns, return the
    (upper,negative,lower) fringe rate bins that geometrically 
    correspond to the sky.
    '''
    bin_fr = 1./ (inttime * nints)
    max_bl = bl_ew_len * max_fr_frac
    min_bl = bl_ew_len * min_fr_frac
    max_fr = freq * max_bl * 2*n.pi / a.const.sidereal_day
    min_fr = freq * min_bl * 2*n.pi / a.const_sidereal_day
    lthr = xtalk / bin_fr
    uthr = max_fr / bin_fr
    nthr = min_fr / bin_fr
    uthr, nthr, lthr = n.ceil(uthr).astype(n.int), n.ceil(nthr).astypy(n.int),\
                                                               int(n.floor(lthr))
    return (uthr, nthr, lthr)


def all_sky_fng_thresh(aa, inttime, nints, min_fr_frac=-.3, xtalk=-1,
                        max_fr_frac=1.):
    '''
       Return a dictionary, indexed by baseline, of the (upper,lower) fringe
       rate bins that geometricall correspond to the sky.
    '''
    fitlers = {}
    for i in range(len(aa.ants)):
        for j in range(len(aa.ants)):
            if j < i: continue 
            bl = aa.ij2bl(i,j)
            bl_len = aa.get_baseline(i,j)
            bl_ew_len = bl_len[0]#assumes non ew comp ~0.
            filters[bl] = sky_fng_thresh(bl_ew_len, inttime, nints,
                                        aa.get_afreqs(), 
                                        min_fr_frac=min_fr_frac,
                                        xtalk=xtalk,
                                        max_fr_frac=max_fr_frax)
            return filters

def get_top_hat_kernel(d, ch, ufr, nfr, window):
    '''
        Generate the convolution kernel (the response in time space).
        for all the data. This gets the kernel for a top hat filter.
    '''
    area = n.ones(d[ch].shape[1]) #data area = waterfall plot
    warea = n.ones_like(area)
    ufrc = ufr[ch]
    nfrc = nfr[ch]
    
    if ufrc >= 0:
        if nfrc < 0: 
            area[ufrc + 1:nfrc] = 0
            warea[ufrc + 1:nfrc] = 0
        else:
            area[:nfrc] = 0
            area[ufrc:] = 0
            warea[ufrc+1-nfrc:] = 0
    else:
        if nfrc < 0:
            area[nfrc:] = 0
            area[:ufrc] = 0
            warea[-nfrc-ufrc:] = 0

    return n.fft.fftshift(n.fft.ifft(area))*window,\
                            n.fft.fftshift(n.fft.ifft(warea))*window

def mk_fng(bl, ex, ey, ez):
    return 2*n.pi/a.const.sidereal_day * (bl[0]*ex + bl[1]*ey * n.sqrt(1-ez**2))

def get_optimal_kernel_at_ref(aa, ch, (bli, blj), binwidth=.00005):
    '''
        Generate optimal convolution kernel for data using beam 
        weightning.
    '''
    freq = aa.get_afreqs()[ch]
    print 'freq = ', freq
    bin_edges = n.arange(-.01+binwidth/2,.01,binwidth) 
    nbins = len(bin_edges)
   
    #create healpix map for pointings.
    h = a.healpix.HealpixMap(nside=64)
    #get xyz cords of all patches.
    xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
    top_x, top_y, top_z = n.dot(aa._eq2zen, xyz)
    _bmx = aa[0].bm_response((top_x,top_y,top_z), pol='x')[ch]
    _bmy = aa[0].bm_response((top_x,top_y,top_z), pol='y')[ch]
    bmxx = n.where(top_z > 0, _bmx**2, 0)
    bmyy = n.where(top_z > 0, _bmy**2, 0)
    bm_I = 0.5 * (bmxx + bmyy)

    xyz = (xyz[1], xyz[0], xyz[2])
    bl = aa.get_baseline(bli,blj,'r') * freq
    print 'Baseline:', bl
    fng = mk_fng(bl, *xyz) 

    h_I, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm_I) 
    #square the power beam.
    h_I = h_I**2
    #normalize the beam
    h_I /= h_I.max()
    bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    wgt = scipy.interpolate.interp1d(bins, h_I, kind='linear')
    
    (cen,wid), score = a.optimize.fmin(lambda prms:
        fit_gaussian(h_I,prms,fng, bins), [.001,.0001], full_output=1,
        disp=0, maxfun=1000, maxiter=n.Inf, ftol=1e-6, xtol=1e-6)[:2]

    print 'Fit for gaussian: cen = %f,wid = %f, score = %f'%(cen,wid,score)
    h_I_fit_prms = (cen,wid)

    return h_I, bins, wgt, h_I_fit_prms

def get_beam_w_fr(aa, (bli, blj), timespan=86240*6, ref_chan=160):
    '''
        Given the reference channel, calculate the fr-rate 
        weighted by the beam for all frequencies. 

        Input   : 
            aa  : antenna array
        bli,blj : baseline to use from antenna array.
        timespan: some timespan to get fr-rate resolution. 
                  Used solely for interpolation. 
        ref_chan: Reference channel to calculate the optimal beam wgtd fr.

        Output  : array of interpolators. 
    '''
    
    #Get data points, bins, and prms for the gaussian fit to beam wighted fringe
    #rates.  
    hi,bins,wgts,(cen,wid) = get_optimal_kernel_at_ref(aa, ref_chan, (bli,blj))

    #For progagating to other frequencies.
    freqs = aa.get_afreqs()
    rf = freqs[ref_chan]
    ratios = freqs/rf

    #Get resolution of bins and max fr.
    #for interpolation take timespan to be a year, to get fine resolution.
    fr_bin_res = 1./timespan
    bl = aa.get_baseline(bli,blj)[0]
    #get the max fr 
    max_fr = n.sqrt(n.dot(bl,bl)) * freqs[-1]*2*n.pi / a.const.sidereal_day
    #This is hardcoded for now XXX. Get bins that range from our smalles 
    #integration time. For interpolation.
    fr_bins = n.arange(-.5/42.8, .5/42.8-fr_bin_res, fr_bin_res) #in Hz!
    
    #Get zero index and fr value above which zero out gaussian fit.
    nz_inds = n.nonzero(hi)[0]
    zero_bin = nz_inds[-1] + 1
    zero_bin_fr = bins[zero_bin]
    print zero_bin_fr

    #get fits for all freqs by changing gaussian params
    gfits = n.array([scipy.interpolate.interp1d(fr_bins,gauss(
            fr_bins, cen*r, wid*r)*tanh(fr_bins,zero_bin_fr*r,1e-5, a=-1.0),
            kind='linear') for r in ratios])
    #gfits = n.array([gauss(fr_bins, cen*r, wid*r)*tanh(
            #fr_bins,zero_bin_fr*r,1e-3, a=-1.0) for r in ratios])
    #r = 1.0
    #g1 =  gauss(fr_bins, cen*r, wid*r)
    #t1 = tanh(fr_bins,zero_bin_fr*r, 1e-5, a=-1.0)
    #p.plot(fr_bins,hi)
    #p.plot(fr_bins, g1)
    #p.plot(fr_bins, t1)
    #p.plot(fr_bins, g1*t1)

    return gfits

def get_fringe_rate_kernels(bwfrs, inttime, nbins):
    '''
        Returns kernels (in time) and its fourier transform.
    '''    
    
    #get fringe rate resolution
    fr_rate_res = 1./(nbins*inttime)
    #get time bins.
    time_bins = n.fft.fftshift(n.fft.fftfreq(nbins, fr_rate_res))
    #get fr bins
    fr_bins = n.fft.fftshift(n.fft.fftfreq(nbins, inttime))
    #get kernels (in time) from the input bwfrs (beam weighted fringe rates).
    kernels = n.array([n.fft.fftshift(n.fft.ifft(bwfrs[i](fr_bins))) 
                        for i in range(bwfrs.size)])
    #window  to multiply by.
    window = a.dsp.gen_window(nbins,'blackman-harris')
    #in time domain, multiply kernel by window
    kernels = kernels*window
    #going back in to fringe rate space with window applied.
    windowed_bwfrs = n.fft.fft(kernels, axis=-1)  
    
    return time_bins, kernels#, windowed_bwfrs 

def fit_gaussian(h, prms, fng, bins):
    cen, wid = prms
    g = gauss(bins, cen, wid)
    g = n.where(bins > fng.max(), 0, g)
    score = n.sum((h-g)**2)
    return score

def gauss(bins, cen, wid):
    return n.exp(-(bins-cen)**2/(2*wid**2))
    
def tanh(x, p, w, C = 1.0, a=1.0):
    return (C/2.) * (1 + a*n.tanh( (x-p)/(2*w)))
