#! /usr/bin/env python
"""
plot the output of SED_fit for a single input trace
"""
import matplotlib
#matplotlib.use('Agg')
import numpy as n,os
import optparse, sys
from pylab import *
import ipdb
from capo.dcj import *
from matplotlib.ticker import NullFormatter
matplotlib.rcParams.update({'font.size':18})

confidence=76.0
gridsize=50
o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for chain in args:
    if not os.path.exists(chain):
        print chain,"not found, skipping"
        continue
    srcname = chain.split('_')[0]
    trace = n.load(chain)
    #compute the parameter confidence limits
    catS0 = find_percentile_errors(10**trace['catlogS0'],confidence)
    catalpha = find_percentile_errors(trace['catalpha'],confidence)
    S0 = find_percentile_errors(10**trace['logS0'],confidence)
    alpha = find_percentile_errors(trace['alpha'],confidence)

    #plot the chains
    figure(1)
    clf()
    subplot(211) #plot the flux chains
    plot(10**trace['catlogS0'],'0.5',alpha=0.5)
    plot(10**trace['logS0'],'k')

    ylabel('flux [Jy]')
    xticks([])
    subplot(212) #plot the SI chains
    plot(trace['catalpha'],'0.5',alpha=0.5)
    plot(trace['alpha'],'k')

    ylabel('Spectral Index')
    xlabel('MCMC step')
    subplots_adjust(wspace=0,left=0.2)
    outfile = srcname + '_traceplot.png'
    print outfile
    savefig(outfile)

    #plot the contours
    figure(2)
    clf()
    left,width = 0.12,0.63
    bottom,height = 0.12,0.63
    bottom_h=left_h = left+width

    rect_contour = [left,bottom,width,height]
    rect_histx = [left,bottom_h,width,0.2]
    rect_histy = [left_h,bottom,0.2,height]
    axContour = axes(rect_contour)
    axHistx = axes(rect_histx)
    axHisty = axes(rect_histy)
   
    if True: #kill the top and left tick labels
        axHistx.xaxis.set_major_formatter(NullFormatter())
        axHisty.yaxis.set_major_formatter(NullFormatter())
    else:
        axHisty.yaxis.tick_right()
        axHistx.xaxis.tick_top()
    if True: #put tick labels for the probability
        axHisty.xaxis.tick_top()
        axHisty.xaxis.set_label_position('top')
        axHistx.yaxis.tick_right()
        axHistx.yaxis.set_label_position('right')
    else:
        axHistx.yaxis.set_major_formatter(NullFormatter())
        axHisty.xaxis.set_major_formatter(NullFormatter())




    #compute the 2D histograms
    pbuf = 1.6
    alpha_lim = (n.min(catalpha+alpha),n.max(catalpha+alpha))
    dalpha = n.max(n.abs([n.diff(catalpha),n.diff(alpha)]))
    S0_lim = (n.min(catS0+S0),n.max(catS0+S0))
    dS0 = n.max(n.abs([n.diff(catS0),n.diff(S0)]))
    alphamin,alphamax  = (alpha_lim[0]-dalpha*pbuf,
        alpha_lim[1]+dalpha*pbuf)
    alphabins = n.linspace(alphamin,
        alphamax,num=gridsize)
    Smin = n.max([S0_lim[0]-dS0*pbuf,0])#don't plot S0 less than 0!
    Smax = S0_lim[1]+dS0*pbuf
    S0bins = n.linspace(Smin,Smax,num=gridsize)
    Hcat,alphabins_cat,S0bins_cat = np.histogram2d(trace['catalpha'],10**(trace['catlogS0']),
                    bins=[alphabins,S0bins])
    H,alphabins,S0bins = np.histogram2d(trace['alpha'],10**(trace['logS0']),
                    bins=[alphabins,S0bins])    
    Pcat = Hcat.astype(n.float)/Hcat.max()
    P = H.astype(n.float)/H.max()
    Plevels = (1 - n.arange(0.1,0.9,0.2))
    Cs = axContour.contour(alphabins[1:],S0bins[1:],P.T,Plevels,colors='k',lw=4)
    axContour.contour(alphabins_cat[1:],S0bins_cat[1:],Pcat.T,Plevels,colors='grey')
    axContour.set_xlim([alphamin,alphamax])
    axContour.set_ylim([Smin,Smax])

    #plot the histograms
    Shist,Sbins = n.histogram(10**trace['logS0'],bins=S0bins)
    Shist_cat,Sbins_cat = n.histogram(10**trace['catlogS0'],bins=S0bins)
    axHisty.plot(Shist.astype(n.float)/Shist.max(),Sbins[1:],'k',lw=4)
    axHisty.plot(Shist_cat.astype(n.float)/Shist_cat.max(),Sbins_cat[1:],'0.5')

    SIhist,SIbins = n.histogram(trace['alpha'],bins=alphabins)
    SIhist_cat,SIbins_cat = n.histogram(trace['catalpha'],bins=alphabins)
    axHistx.plot(SIbins[1:],SIhist.astype(n.float)/SIhist.max(),'k',lw=4)
    axHistx.plot(SIbins_cat[1:],SIhist_cat.astype(n.float)/SIhist_cat.max(),'0.5')

    #plot the 76% confidence lines

    
    axHisty.set_ylim([Smin,Smax])
    axHisty.set_xlim([0.01,0.99])
    axHisty.set_xlabel('$\mathtt{P}(M|D)$',size=10)
    axHistx.set_ylabel('$\mathtt{P}(M|D)$',size=10)
    setp(axHisty.get_xticklabels(),fontsize=12)
   

    axHistx.set_xlim([alphamin,alphamax])
    axHistx.set_ylim([0.01,0.99])
    setp(axHistx.get_yticklabels(),fontsize=12)

    axContour.set_ylabel('$S150$ [Jy]')
    axContour.set_xlabel('$\\alpha$')


    outfile = srcname+'_trace_hist.png'
    print outfile
    savefig(outfile)
