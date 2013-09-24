#! /usr/bin/env python
"""
Produce mcmc contours comparing the calibration fit contour with 
the regular specfind fits.
Run after fit_calibrators and SED_fit

usage
plot_cal_contours.py srcA_srcB_gain_chain.npz
"""
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from capo.dcj import *
from pylab import *
from scipy import optimize
import ipdb
matplotlib.rcParams.update({'font.size':18})
confidence=73.
pbuf = 1.6
gridsize=50
caltrace = sys.argv[-1]
cals = caltrace.split('_gain')[0].split('_')
print "loading %d sources"%(len(cals))
D = n.load(caltrace)['flatchain']
for i,srcname in enumerate(cals):
    #load the cal trace
    S_gtrace = 10**D[:,2*i+1]
    alpha_gtrace = D[:,2*i]
    #compute the confidence intervals
    gS0 = find_percentile_errors(S_gtrace,confidence) #with more cal points
    galpha = find_percentile_errors(alpha_gtrace,confidence)
    print a2l(gS0)
    print a2l(galpha)

    #load the ned chain 
    chainout = srcname+'_mcmc_chain.npz'
    trace = n.load(chainout)
    if not os.path.exists(chainout):
#        print "mcmc trace file %s not found, skipping analysis step"%(chainout+'.npz')
        continue
    S0 = find_percentile_errors(10**trace['logS0'],confidence) #with specfind + PAPER
    alpha = find_percentile_errors(trace['alpha'],confidence)
    catS0 = find_percentile_errors(10**trace['catlogS0'],confidence) #with just specfind
    catalpha = find_percentile_errors(trace['catalpha'],confidence)
    
    #load the specfind chain
    chainout_sf = srcname+'_mcmc_chain_specfind.npz'
    trace_sf = n.load(chainout_sf)
    S0_sf = find_percentile_errors(10**trace_sf['logS0'],confidence) #with specfind + PAPER
    alpha_sf = find_percentile_errors(trace_sf['alpha'],confidence)
    catS0_sf = find_percentile_errors(10**trace_sf['catlogS0'],confidence) #with just specfind
    catalpha_sf = find_percentile_errors(trace_sf['catalpha'],confidence)

    #decide what bin ranges to use based on all 5 traces
    alpha_lim = (n.min(catalpha+alpha+galpha+alpha_sf+catalpha_sf),
        n.max(catalpha+alpha+galpha+alpha_sf+catalpha_sf))
    dalpha = n.max(n.abs([n.diff(catalpha),n.diff(alpha)]))
    S0_lim = (n.min(catS0+S0+gS0+S0_sf+catS0_sf),n.max(catS0+S0+gS0+S0_sf+catS0_sf))
    dS0 = n.max(n.abs([n.diff(catS0),n.diff(S0)]))
    alphabins = n.linspace(alpha_lim[0]-dalpha*pbuf,
        alpha_lim[1]+dalpha*pbuf,num=gridsize)
    S0bins = n.linspace(S0_lim[0]-dS0*pbuf,S0_lim[1]+dS0*pbuf,num=gridsize)

    #compute the histograms
    #the ned fit histograms
    Hcat,alphabins_cat,S0bins_cat = np.histogram2d(trace['catalpha'],10**(trace['catlogS0']),
            bins=[alphabins,S0bins])
    Pcat = Hcat/Hcat.max()
    H,alphabins,S0bins = np.histogram2d(trace['alpha'],10**(trace['logS0']),
            bins=[alphabins,S0bins])
    P = H/H.max()
    #the gain fit histograms
    gH,galphabins,gS0bins = np.histogram2d(alpha_gtrace,S_gtrace,
            bins=[alphabins,S0bins])
    gP = gH/gH.max()
    #the specfind histograms
    Hcat_sf,alphabins_cat_sf,S0bins_cat_sf = np.histogram2d(trace_sf['catalpha'],10**(trace_sf['catlogS0']),
            bins=[alphabins,S0bins])
    Pcat_sf = Hcat_sf/Hcat_sf.max()
    H_sf,alphabins_sf,S0bins_sf = np.histogram2d(trace_sf['alpha'],10**(trace_sf['logS0']),
            bins=[alphabins,S0bins])
    P_sf = H_sf/H_sf.max()

    figure(1)
    clf()
    title(srcname)
    #plot the specfind contours
    contour(alphabins[1:],S0bins[1:],P_sf.T,[P_sf.max()*(1-confidence/(100.))],
                colors='k')
    contour(alphabins_cat[1:],S0bins_cat[1:],Pcat_sf.T,[
                                Pcat_sf.max()*(1 - confidence/(100.)),
                                                    ],colors='grey')
    xlabel('Spectral Index')
    ylabel('$S_{150}$ [Jy]')
    savefig(srcname+'_SI_MCMC_specfind.png')
    #plot the ned contours
    contour(alphabins[1:],S0bins[1:],P.T,[P.max()*(1-confidence/(100.))],
                colors='k')
    contour(alphabins_cat[1:],S0bins_cat[1:],Pcat.T,[                                                                Pcat.max()*(1 - confidence/(100.)),
                                                    ],colors='grey')
    savefig(srcname+'_SI_MCMC_specfind_ned.png')
    ylim([S0_lim[0]-dS0*pbuf,S0_lim[1]+dS0*pbuf]) #awful machinations to get axis limits to come out nice
    xlim([alpha_lim[0]-dalpha*pbuf,
            alpha_lim[1]+dalpha*pbuf]) #but not backwards sometimes      
    contour(alphabins[1:],S0bins[1:],gP.T,[gP.max()*(1 - confidence/100.)],
                    colors='k',linestyles='dotted')

    savefig(srcname+'_SI_MCMC_specfind_ned_cal.png')


    
