#! /usr/bin/env python
"""
Plot the ouput of beam_src_vs_ha along with other catalog data
"""
import matplotlib,emcee
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from capo.dcj import *
from pylab import *
from scipy import optimize
import ipdb
matplotlib.rcParams.update({'font.size':18})
CAT=True
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
o.add_option('--mcmc_recompute',action='store_true',
    help='Ignore saved traces and rerun mcmc fit, otherwise just move on to plotting.')
o.add_option('--plot',action='store_true',
    help='Do MCMC contour plots')
o.add_option('--calmodel',
    help="76% confidence output of catalog_spectrum_fit for calibrator.format=srcname S0_low,S0_val,S0_high alpha_low,alpha_val,alpha_high ")
o.add_option('-v',action='store_true',
    help='output more stuff')
o.add_option('--fmax',default=5e3,type='float',
    help='Maximum catalog frequency in MHz [default=5e3]')
o.add_option('--nsteps',default=1000,type=n.int,
    help='Number of MCMC steps [default=1000]')
opts,args = o.parse_args(sys.argv[1:])


#Set some constants 
confidence = 73.
#set up the chain
nwalkers = 100
nsteps = opts.nsteps
model_select='SI'

#some hardcoded plotting options
rows=6
cols=4
PLOT_CONFIDENCE=opts.plot
PLOT_CONFIDENCE_subplot=False
alpha_range = [-2,0]
S0rangef = n.array([.95,1.05]) #this in percent of median flux
gridsize = 50 #number of bins in probability histogram
#define the likelihood model
#load the specfind catalog
#load the new points
#for each new point, fit a model and output a trace
#if a trace is found and we aren't recomputing, skip the mcmc and analyze the trace
#analysis should be peak location (S,SI), error 76% confidence interval on same


"""
Define the probability
So far I am assuming the flux error bars are 1D and we can use the usual gaussian likelihood"""
def lnSI_lin(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    S0 = 10**x[1]
    return  -1 * n.sum( ( fluxes - S0 * (freqs/f0)**alpha)**2 / (2 * errors**2) )
def lnSI_exp(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = (logS_high-logS_low)/2#the error in log space
    logS = n.log10(fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0)**2/(2*dlogS**2))#the log likelihood
def lnSI_exp_wgain(x,freqs,fluxes,errors,gainchain,ispaper):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = (logS_high-logS_low)/2#the error in log space
    G = n.ones_like(fluxes)
    G[ispaper>0] *= n.random.choice(gainchain)
    logS = n.log10(G*fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0)**2/(2*dlogS**2))#the log likelihood

def lnSI_wgain(x,freqs,fluxes,errors,gainchain,ispaper):
    f0=150.
    alpha = x[0]
    S0 = 10**x[1]
    G = n.ones_like(fluxes)
    G[ispaper>0] *= n.random.choice(gainchain)
    return  -1 * n.sum( ( G*fluxes - S0 * (freqs/f0)**alpha)**2 / (2 * errors**2) )
    
lnSI = lnSI_exp
def gain_errors(calmodel,freqs):
    #parse the calibrator model
    S0_range = map(float,calmodel.split()[1].split(','))
    alpha_range = map(float,calmodel.split()[2].split(','))
    dS0 = (S0_range[-1] - S0_range[0])/2
    S0 = S0_range[1]
    dalpha = (alpha_range[-1] - alpha_range[0])/2
    alpha = alpha_range[1]
    return n.abs(dS0)/S0 * n.ones_like(freqs) + n.abs(alpha * dalpha/(freqs/150.))


def pic_spectrum(nedxmlfile):
    ############
    ### set some data from NED
    ####
    NED = atpy.Table('pic_ned_spectrum.vot')
    fluxcol = 'photo_col7'
    freqcol = 'photo_col6'
    errorcol = 'photo_col8'
    units = 'photo_col5'
    NEDfreqs = NED[freqcol]
    NED = NED.where(NEDfreqs<5000e6)
    nedfreqs = n.array(NED[freqcol]/1e6)
    nedfluxes = n.array(NED[fluxcol])
    nederrorstrings = NED[errorcol]
    nederrors = n.array([0]*len(nedfluxes))
    for i,err in enumerate(nederrorstrings):
        if len(err)<2: continue
        try:
            nederrors[i] = float(err.split('+/-')[1])
        except(ValueError): 
            if err.split('+/-')[1].endswith('%'):
                nederrors[i] = nedfluxes[i]*float(err.split('+/-')[1].split('%')[0].strip())/100
    return nedfreqs[nederrors>0],nedfluxes[nederrors>0],nederrors[nederrors>0]
def load_paper_spectra(filename):
    lines = open(filename).readlines()
    spectra ={}
    for line in lines:
        if line.startswith('#'):
            freqs = n.array(map(float,line.split('=')[1].split(','))).squeeze()
            continue
        line = line.split()
        srcname = line[0]
        fluxes = n.array(map(float,line[1].split(',')))
        errors = n.array(map(float,line[2].split(',')))
        spectra[srcname] = {'flux':fluxes,'err':errors}
    return freqs,spectra



#load teh new points
PAPERfreqs,PAPERcatalog = load_paper_spectra(args[0])
dG = n.zeros_like(PAPERfreqs)
gainchain = n.ones(10)
G = 1.0
cals = []
#check for the calibrator and calculate the gain errors
if not opts.calmodel is None:
    if not os.path.exists(opts.calmodel):
    #   errors due to uncertainty in flux calibrator catalog values
    #   calculated by fitting with MCMC method using catalog_spectrum_fit.py
        dG = gain_errors(opts.calmodel,PAPERfreqs)
        print "median catalog calibrator uncertainty gain error %4.1f%%"%(n.median(dG)*100)
        print "min/max %4.1f%%/%4.1f%%"%(dG.min()*100,dG.max()*100)
        try:
            dG += PAPERcatalog[opts.calmodel.split()[0]]['err']/PAPERcatalog[opts.calmodel.split()[0]]['flux']
        except(KeyError):
            print "WARNING"
            print "calibrator %s not found in input catalog %s. Calibration errors cannot be factored in!!"
        print "#median gain error %4.1f%%"%(n.median(dG)*100)
        print "#min/max gain error  %4.1f%%/%4.1f%%"%(dG.min()*100,dG.max()*100)
    elif os.path.exists(opts.calmodel):
        #get a list of the calibrators
        try:
            cals = opts.calmodel.split('gain')[0].split('_')
        except(IndexError):
            pass
        #load a posterior sample of the gains
        print "#using gain tracefile ",opts.calmodel
        gainchain = n.load(opts.calmodel)['flatchain'][:,-1]
        gain = n.array(find_percentile_errors(dB(n.abs(gainchain)),confidence))
        print "#gain [dB]", a2l(n.round(gain,2))
        print "#gain [x]",a2l(n.round(idB(gain),2))
        gain = idB(gain)
        dG += (gain[-1] - gain[0])/2/gain[1]
        G = gain[1]       
        print "#gain error [mult,%]", (gain[-1] - gain[0])/2/gain[1]*100    
    print   "GAIN MODEL:",G,'+/-',dG
#load the known spectral data
if CAT:
    """
    Loading the whole catalog takes forever. If a cache exists. USE IT!
    """
    specfind_cache_file = 'specfind_'+args[0][:-4]+'.vot'
    if not os.path.exists(specfind_cache_file):
        print 'specfind selection cache ',specfind_cache_file,'not found. Generating...'
        specfind = atpy.Table('specfind_spectra.vot')
        for i,srcname in enumerate(PAPERcatalog):
            if i==0: specfind_cache = select_table_where_source(specfind,srcname)
            else:
                selection = select_table_where_source(specfind,srcname)
                specfind_cache.append(selection)
        specfind_cache.write(specfind_cache_file)
    else:
        print "loading specfind subset cache: ",specfind_cache_file
        specfind = load_table(specfind_cache_file)
        
            
PAPERcatalog.pop('0518-458B')
#for each source, do the fit
badfits = []
for srcname in sort(PAPERcatalog.keys()):
    chainout = srcname+'_mcmc_chain.npz'
    if srcname =='0518-458B':continue #this is just pictor all over again. skip
    if os.path.exists(chainout) and not opts.mcmc_recompute:
#        print "mcmc trace file %s exists, skipping"%(chainout)
        continue
    print "MCMC on ",srcname

    fluxes,errors = PAPERcatalog[srcname]['flux'],PAPERcatalog[srcname]['err']
    fluxes *= G
    errors += PAPERcatalog[srcname]['flux'] * dG #factor in the calibration error
    #remove negative data points
    pos = n.argwhere(fluxes>0).squeeze()
    freqs,fluxes,errors = PAPERfreqs[pos],fluxes[pos],errors[pos]
    if CAT:
        nedfile = srcname+'_ned_spectrum.vot'
        if os.path.exists(nedfile):
            print "plotting ned data in:",nedfile
            cat_freq,cat_flux,cat_err = ned_spectrum(nedfile,doprint=opts.v,fmax=opts.fmax*1e6)
        else:
            print "using specfind catalog"
            cat_table = select_table_where_source(specfind,srcname)
            cat_table.sort('nu')
            cat_freq,cat_flux,cat_err = spectrum(cat_table,fmax=opts.fmax)
        #fold in the catalog data
        try:            
            if cat_flux.size>0:
                if cat_flux.size==1:
                    cat_freq.shape,cat_flux.shape,cat_err.shape = (1,)*3
                else:
                    cat_freq = cat_freq.squeeze()
                print "folding in %d data points from specfind catalog"%cat_flux.size 
                print cat_freq,cat_flux,cat_err
                freqs =  n.concatenate((freqs,cat_freq))
                fluxes = n.concatenate((fluxes,cat_flux)) 
                errors = n.concatenate((errors,cat_err))
                ispaper = n.concatenate((n.ones_like(freqs),n.zeros_like(cat_err)))
            else:
                print "not enough data found in specfind table"
        except(TypeError):
            print "not enough data found in specfind table"
            continue
        print freqs,fluxes,errors
    """
    FITTING THE MODEL
    """
    #ipdb.set_trace() 
    if model_select == 'SI': 
        print "#using gaussian model"
        ndim = 2
        #set up the initial positions
        p0 = np.vstack([ #regular 2d error model
            np.random.uniform(-1,1,size=nwalkers), #alpha
            np.random.uniform(1,2,size=nwalkers),   #S0
            ]).T
        if not opts.calmodel is None and os.path.exists(opts.calmodel):
            sampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,
#                args=[freqs,fluxes,errors,gainchain,ispaper])
                args=[freqs,fluxes,errors])
        else:
            sampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,
                args=[freqs,fluxes,errors])  
        if cat_flux.size>1:
            catsampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,args=[cat_freq,cat_flux,cat_err]) 
    
    print "#burning in"
    #run the sampler 100 steps for burn-in
    pos, prob, state = sampler.run_mcmc(p0, 100)
    if cat_flux.size>1:
        pos, prob, state = catsampler.run_mcmc(p0, 100)
    
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    if cat_flux.size>1:
        catsampler.reset()

    # Starting from the final position in the burn-in chain, sample for nsteps
    # steps.
    print "#sampling witn %d walkers for a %d steps"%(nwalkers,nsteps)
    if cat_flux.size>1:
        catsampler.run_mcmc(pos,nsteps,rstate0=state)   
        print "Mean acceptance fraction (Catalog):",np.mean(catsampler.acceptance_fraction)

    sampler.run_mcmc(pos,nsteps,rstate0=state)
    print "#Mean acceptance fraction (PAPER):", np.mean(sampler.acceptance_fraction)
    if np.mean(sampler.acceptance_fraction)<0.25: badfits.append(srcname)
    alpha = sampler.flatchain[:,0]
    logS0 = sampler.flatchain[:,1]
    if n.median(logS0)>100: ipdb.set_trace()
    if cat_flux.size>1:
        np.savez(chainout,alpha=alpha,logS0=logS0,flatchain=sampler.flatchain,catalpha=catsampler.flatchain[:,0],catlogS0=catsampler.flatchain[:,1])
    else:
        np.savez(chainout,alpha=alpha,logS0=logS0,flatchain=sampler.flatchain,catalpha=n.zeros(nsteps*nwalkers),catlogS0=n.zeros(nsteps*nwalkers))
#for each source do the histogram of the trace.
print "#Name Flux_limits SI_limits Catalog_Flux_limits Catalog_SI_limits"
for i,srcname in enumerate(sort(PAPERcatalog.keys())):
    if srcname in ['0518-458B']:continue #this is just pictor all over again. skip
#    print "MCMC on ",srcname
    chainout = srcname+'_mcmc_chain.npz'
    if not os.path.exists(chainout):
#        print "mcmc trace file %s not found, skipping analysis step"%(chainout+'.npz')
        continue
    trace = n.load(chainout)
    goodcatfit = n.sum(trace['catlogS0'])>0
    n.set_printoptions(precision=2)
    #print "\tS0",
    S0 = find_percentile_errors(10**trace['logS0'],confidence)
    alpha = find_percentile_errors(trace['alpha'],confidence)
    

    if goodcatfit: 
        catS0 = find_percentile_errors(10**trace['catlogS0'],confidence)
        catalpha = find_percentile_errors(trace['catalpha'],confidence)
    if n.max(find_percentile_errors(trace['logS0'],confidence))>100:
#        print 'BAD!'
        continue
    else:
        if srcname in badfits: print '#'+srcname
        else:print srcname,
        #print the PAPER+CAT flux conf
        print a2l(n.round(S0,2)),
#        print "\talpha",
        print a2l(n.round(find_percentile_errors(trace['alpha'],confidence),2)),
#        print "\t alpha catalog",
        #print the CAT flux conf
        if goodcatfit:
            pbuf = 1.6
            alpha_lim = (n.min(catalpha+alpha),n.max(catalpha+alpha))
            dalpha = n.max(n.abs([n.diff(catalpha),n.diff(alpha)]))
            S0_lim = (n.min(catS0+S0),n.max(catS0+S0))
            dS0 = n.max(n.abs([n.diff(catS0),n.diff(S0)]))
            alphabins = n.linspace(alpha_lim[0]-dalpha*pbuf,
                alpha_lim[1]+dalpha*pbuf,num=gridsize)
            S0bins = n.linspace(S0_lim[0]-dS0*pbuf,S0_lim[1]+dS0*pbuf,num=gridsize)

            print a2l(n.round(catS0,2)),
            print a2l(n.round(find_percentile_errors(trace['catalpha'],confidence),2)),
            Hcat,alphabins_cat,S0bins_cat = np.histogram2d(trace['catalpha'],10**(trace['catlogS0']),
                    bins=[alphabins,S0bins])
            H,alphabins,S0bins = np.histogram2d(trace['alpha'],10**(trace['logS0']),
                    bins=[alphabins,S0bins])
        else:
            print "-1,-1,-1",
            print "-99,-99,-99",
            H,alphabins,S0bins = np.histogram2d(trace['alpha'],10**(trace['logS0']),
                    bins=gridsize)

        #print the FOM for PAPER+CAT
        P = H/H.sum()
        FOM = 1./n.sum(P>(P.max()*(1-confidence/100)))*100
        print n.round(FOM,2),
        #print the FOM for CAT
        if goodcatfit:
            Pcat = Hcat/Hcat.sum()
            FOMcat = 1./n.sum(Pcat>(Pcat.max()*(1-confidence/100)))*100
            print n.round(FOMcat,2),
        #print the improvement in FOM
            print n.round(FOM-FOMcat,2),
            P76 = n.round(P/(P.max()*(1-confidence/100)))
            Pcat76 = n.round(Pcat/(Pcat.max()*(1-confidence/100)))
            overlap = n.sum(n.logical_and(P76,Pcat76))/n.sum(P76) #fractional contour overlap
            catarea = n.sum(Pcat76)

            #print n.round(n.dot(Pcat.ravel(),P.ravel())*1e4,2),#/(n.sum(P)*n.sum(Pcat))
            print n.round(overlap,2),
            print n.round((FOM-FOMcat)*overlap,4),
            print n.sum(n.abs(n.diff(P76)))
            #for large improvements to existing models overlap<<catarea
        else:
            print 0,FOM,0,0,n.sum(n.abs(n.diff(n.round(P/(P.max()*(1-confidence/100))))))
    if PLOT_CONFIDENCE:
        print "#plotting confidence"

        figure(2)
        clf()
        #print "#extent = ",n.round(extent,4)
        #print "#PMAX = ",P.max()
        #print "alpha_lim,S0_lim",n.round(alpha_lim,4),n.round(S0_lim,4)
        contour(alphabins[1:],S0bins[1:],P.T,[#P.max()*99.99/(2*100.),
                                                                P.max()*(1-confidence/(100.)),
                                                               # P.max()*50/(100.)
                                                                ],colors='k')
        #compute prior catalog histogram
        if goodcatfit:
            contour(alphabins_cat[1:],S0bins_cat[1:],Pcat.T,[#P.max()*99.99/(2*100.),
                                                                Pcat.max()*(1 - confidence/(100.)),
                                                                #Pcat.max()*50/(2*100.)
                                                                ],colors='grey')

            ylim([S0_lim[0]-dS0*pbuf,S0_lim[1]+dS0*pbuf]) #aweful machinations to get axis limits to come out nice
            xlim([alpha_lim[0]-dalpha*pbuf,
            alpha_lim[1]+dalpha*pbuf]) #but not backwards sometimes                                                              
        #annotate(srcname,[0.5,0.85],xycoords='axes fraction',textcoords='axes fraction',size=12)
        if srcname in cals: title(srcname+'*')
        else:title(srcname)
        #plot the previous known values
#        cat_table = select_table_where_source(specfind,srcname)
#        plot(cat_table['a'],cat_table['S_nu_']/1e3*(150./cat_table['nu'])**cat_table['a'],'s')
        xlabel('Spectral Index')
        ylabel('$S_{150}$ [Jy]')
        #savefig(srcname+'_SI_MCMC.eps')
        savefig(srcname+'_SI_MCMC.png')


            
print "BADFITS"
print "#\n".join(badfits)
                


