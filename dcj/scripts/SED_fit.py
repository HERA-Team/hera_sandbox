#! /usr/bin/env python
"""
Plot the ouput of beam_src_vs_ha along with other catalog data
"""
import matplotlib,emcee
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from pylab import *
from scipy import optimize
import ipdb

CAT=True
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
o.add_option('--mcmc_recompute',action='store_true',
    help='Ignore saved traces and rerun mcmc fit, otherwise just move on to plotting.')
o.add_option('--plot',action='store_true',
    help='Do MCMC contour plots')
opts,args = o.parse_args(sys.argv[1:])


#Set some constants 
confidence = 73.
#set up the chain
nwalkers = 100
nsteps = 1000
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
def lnSI(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = logS_high-logS_low#the error in log space
    logS = n.log10(fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0)**2/(2*dlogS**2))#the log likelihood
def select_table_where_source(table,srcname):
    #Assumes MRC names!!!
    #make a special exception for pic
    if srcname=='pic':
        seq=20658
    elif srcname=='cen':
        seq=61767
    else:
        result = table.where(table['Name']== 'MRC %s'%srcname)
        seq = result['Seq']
    return table.where(table['Seq']==seq)
def spectrum(table):
    freqs,fluxes,errors = n.array(table['nu']),n.array(table['S_nu_']/1e3),n.array(table['e_S_nu_']/1e3)
    try: len(freqs)
    except(TypeError):
        freqs.shape = (1,)
        fluxes.shape = (1,)
        errors = (1,)
    return freqs,fluxes,errors
def pic_spectrum(nedxmlfile):
    ############
    ### set some data from NED
    ####
    NED = atpy.Table('pic_ned_spectrum.xml')
    fluxcol = 'photo_col7'
    freqcol = 'photo_col6'
    errorcol = 'photo_col4'
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
def find_trace_peak(trace,gridsize):
    #input trace and number of grid points
    #output indices and values
    H,bins = np.histogramdd(trace)
    peak_index = H.argmax()
    return peak_index,[bins[i][peak] for i,j in enumerate(peak_index)]
def find_closest(A,a):
    return np.abs(A-a).argmin()
def find_percentile_errors(trace,percentile,nbins=100):
    thist,bins = np.histogram(trace,bins=nbins)
    binwidth = np.diff(bins)[0]
    lower = bins[find_closest(np.cumsum(thist)/float(np.sum(thist)),
        0.5-(percentile/100.)/2)]+binwidth/2.
    upper = bins[find_closest(np.cumsum(thist)/float(np.sum(thist)),
        0.5+(percentile/100.)/2)]+binwidth/2.
    med = np.median(trace)
    return [lower,med,upper]
def get_votable_column_names(tablefile,tid=0):
    from atpy.votable import parse
    votable = parse(tablefile)
    colnames = {}
    for id, table in enumerate(votable.iter_tables()):
        if id==tid:
            break    
    for field in table.fields:
        colnames[field._ID] = field.name
    return colnames
def load_table(tablefile):
    colnames = get_votable_column_names(tablefile)
    table = atpy.Table(tablefile)
    for ID,name in colnames.iteritems():
        table.rename_column(ID,name)
    return table
def a2l(array):
    return ','.join(map(str,array))


#load teh new points
PAPERfreqs,PAPERcatalog = load_paper_spectra(args[0])
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
for srcname in sort(PAPERcatalog.keys()):
    chainout = srcname+'_mcmc_chain.npz'
    if srcname =='0518-458B':continue #this is just pictor all over again. skip
    if os.path.exists(chainout) and not opts.mcmc_recompute:
#        print "mcmc trace file %s exists, skipping"%(chainout)
        continue
    print "MCMC on ",srcname

    fluxes,errors = PAPERcatalog[srcname]['flux'],PAPERcatalog[srcname]['err']
    if CAT:
        #load the source spectrum from the prior data catalog
        if srcname!='pic':
            print "using specfind catalog"
            cat_table = select_table_where_source(specfind,srcname)
            cat_freq,cat_flux,cat_err = spectrum(cat_table)
        else:
            print "using NED data"
            cat_freq,cat_flux,cat_err = pic_spectrum('pic_ned_spectrum.xml')

        #print "folding in %d data points from specfind catalog"%cat_flux.size
        
        #fold in the catalog data
        try:
            print "folding in %d data points from specfind catalog"%cat_flux.size
            if cat_flux.size>1:
                cat_freq = cat_freq.squeeze()
            print cat_freq,cat_flux,cat_err
            freqs =  n.concatenate((PAPERfreqs,cat_freq))
            fluxes = n.concatenate((fluxes,cat_flux)) 
            errors = n.concatenate((errors,cat_err))
        except(TypeError):
            print "not enough data found in specfind table"
            continue
        print freqs,fluxes,errors
    """
    FITTING THE MODEL
    """
    #ipdb.set_trace() 
    if model_select == 'SI': 
        print "using gaussian model"
        ndim = 2
        #set up the initial positions
        p0 = np.vstack([ #regular 2d error model
            np.random.uniform(-1,1,size=nwalkers), #alpha
            np.random.uniform(1,2,size=nwalkers),   #S0
            ]).T
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,args=[freqs,fluxes,errors])  
        catsampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,args=[cat_freq,cat_flux,cat_err]) 
    
    print "burning in"
    #run the sampler 100 steps for burn-in
    pos, prob, state = sampler.run_mcmc(p0, 100)
    pos, prob, state = catsampler.run_mcmc(p0, 100)
    
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    catsampler.reset()
    # Starting from the final position in the burn-in chain, sample for nsteps
    # steps.
    print "sampling witn %d walkers for a %d steps"%(nwalkers,nsteps)
    sampler.run_mcmc(pos,nsteps,rstate0=state)
    catsampler.run_mcmc(pos,nsteps,rstate0=state)
    print "Mean acceptance fraction:", np.mean(sampler.acceptance_fraction),np.mean(catsampler.acceptance_fraction)
    
    alpha = sampler.flatchain[:,0]
    logS0 = sampler.flatchain[:,1]
    if n.median(logS0)>100: ipdb.set_trace()
    np.savez(chainout,alpha=alpha,logS0=logS0,flatchain=sampler.flatchain,catalpha=catsampler.flatchain[:,0],catlogS0=catsampler.flatchain[:,1])
#for each source do the histogram of the trace.
print "#Name Flux_limits SI_limits Catalog_Flux_limits Catalog_SI_limits"
for i,srcname in enumerate(sort(PAPERcatalog.keys())):
    if srcname =='0518-458B':continue #this is just pictor all over again. skip
#    print "MCMC on ",srcname
    chainout = srcname+'_mcmc_chain.npz'
    if not os.path.exists(chainout):
#        print "mcmc trace file %s not found, skipping analysis step"%(chainout+'.npz')
        continue
    trace = n.load(chainout)
    n.set_printoptions(precision=2)
    #print "\tS0",
    logS0 = find_percentile_errors(trace['logS0'],confidence)
    catlogS0 = find_percentile_errors(trace['catlogS0'],confidence)
    if n.max(logS0)>100:
#        print 'BAD!'
        continue
    else:
        print srcname,
        #print the PAPER+CAT flux conf
        print a2l(n.round(10**n.array(logS0),2)),
#        print "\talpha",
        print a2l(n.round(find_percentile_errors(trace['alpha'],confidence),2)),
#        print "\t alpha catalog",

        #print the CAT flux conf
        print a2l(n.round(10**n.array(catlogS0),2)),
        print a2l(n.round(find_percentile_errors(trace['catalpha'],confidence),2)),

        #print the FOM for PAPER+CAT
        H,alphabins,S0bins = np.histogram2d(trace['alpha'],10**(trace['logS0']),
                bins=gridsize)
        P = H/H.sum()
        FOM = 1./n.sum(P>(P.max()*confidence/(2*100)))*100
        print n.round(FOM,2),
        #print the FOM for CAT
        Hcat,alphabins,S0bins = np.histogram2d(trace['catalpha'],10**(trace['catlogS0']),
                bins=[alphabins,S0bins])
        Pcat = Hcat/Hcat.sum()
        FOMcat = 1./n.sum(Pcat>(Pcat.max()*confidence/(2*100)))*100
        print n.round(FOMcat,2),
        #print the improvement in FOM
        print n.round(FOM-FOMcat,2)
    if PLOT_CONFIDENCE:
        figure(2)
        clf()
        print "computing histogram"
        contour(alphabins[1:],S0bins[1:],P.T,[#P.max()*99.99/(2*100.),
                                                                P.max()*confidence/(2*100.),
                                                                #P.max()*50/(2*100.)
                                                                ],colors='k')
        #compute prior catalog histogram
        contour(alphabins[1:],S0bins[1:],Pcat.T,[#P.max()*99.99/(2*100.),
                                                                Pcat.max()*confidence/(2*100.),
                                                                #P.max()*50/(2*100.)
                                                                ],colors='grey')


        #ylim(S0range)
        #xlim(alpha_range)                                                              
        annotate(srcname,[0.5,0.75],xycoords='axes fraction',textcoords='axes fraction',size=9)
        #plot the previous known values
        cat_table = select_table_where_source(specfind,srcname)
        plot(cat_table['a'],cat_table['S_nu_']/1e3*(150./cat_table['nu'])**cat_table['a'],'s')
        xlabel('spectral index')
        ylabel('flux @ 150MHz [Jy]')
        savefig(srcname+'_SI_MCMC.eps')
        savefig(srcname+'_SI_MCMC.png')

        if False:
            figure()
            clf()
                    #compute prior catalog histogram
            Hcat,alphabins,S0bins = np.histogram2d(trace['catalpha'],10**(trace['catlogS0']),
                    bins=gridsize)
            Pcat = Hcat/Hcat.sum()
            contour(alphabins[1:],S0bins[1:],Pcat,[#P.max()*99.99/(2*100.),
                                                                    P.max()*confidence/(2*100.),
                                                                    #P.max()*50/(2*100.)
                                                                    ],colors='grey')
            xlabel('spectral index')
            ylabel('flux @ 150MHz [Jy]')
            savefig(srcname+'cat_SI_MCMC.eps')
            savefig(srcname+'cat_SI_MCMC.png')
            

                


