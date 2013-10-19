#! /usr/bin/env python
"""
fit a spectral model to a single source. 
used for fitting calibrator model in psa64 beamform pic strip paper

usage
catalog_spectrum_fit.py specfind_catalog.vot 2331-416

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
o.add_option('-v',action='store_true',
    help="Print more details about fitting process and catalog parsing")
opts,args = o.parse_args(sys.argv[1:])

catalog = args[0]
srcname = args[1]
#Set some constants 
confidence = 73.
#set up the chain
nwalkers = 100
nsteps = 5000
model_select='SI'

gridsize = 50 #number of bins in probability histogram

def lnSI_lin(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    S0 = 10**x[1]
    return  -1 * n.sum( ( fluxes - S0 * (freqs/f0)**alpha)**2 / (2 * errors**2) )
def lnSI_curved(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    curve = x[2]
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = (logS_high-logS_low)/2#the error in log space
    logS = n.log10(fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0 - curve**2 * n.log10(freqs/f0))**2/(2*dlogS**2))#the log likelihood

def lnSI_exp(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = (logS_high-logS_low)/2#the error in log space
    logS = n.log10(fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0)**2/(2*dlogS**2))#the log likelihood
lnSI = lnSI_curved
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
def ned_spectrum(nedxmlfile,doprint=False):
    ############
    ### set some data from NED
    ####
    NED = atpy.Table(nedxmlfile)
    fluxcol = 'photo_col7'
    freqcol = 'photo_col6'
    errorcol = 'photo_col8'
    unitscol = 'photo_col5'
    refcol = 'photo_col10'
    NEDfreqs = NED[freqcol]
    NED = NED.where(NEDfreqs<5000e6)
    nedfreqs = n.array(NED[freqcol]/1e6)
    nedfluxes = n.array(NED[fluxcol])
    nederrorstrings = NED[errorcol]
    nedrefs = NED[refcol]
    #nederrors = n.array([0.]*len(nedfluxes)) %Danny being conservative. No error bar = somehow suspect measurement.
    #    in this past has indicated duplicates data using an older calibration
    nederrors = nedfluxes*0.25 #arp being liberal, assuming 25% error bars for data with no catalog error
    guess = n.ones_like(nedfluxes).astype(n.bool)
    for i,err in enumerate(nederrorstrings):
        if len(err)<2: continue
        try:
            nederrors[i] = float(err.split('+/-')[1])
            guess[i] = False
        except(ValueError): 
            if err.split('+/-')[1].endswith('%'):
                nederrors[i] = nedfluxes[i]*float(err.split('+/-')[1].split('%')[0].strip())/100
                guess[i] = False
    if doprint:
        print "using %s with the following error bars"%nedxmlfile
        print "freq\tflux\tError\tFrac\tGuess?\tCite"
        for i in range(len(nederrors)):
            print "%7.2f\t%7.2f\t%7.2f\t%7.2f\t%r\t%s"%(
                    nedfreqs[nederrors>0][i],
                    nedfluxes[nederrors>0][i],
                    nederrors[nederrors>0][i],
                    nederrors[nederrors>0][i]/nedfluxes[nederrors>0][i],
                    guess[nederrors>0][i],
                    nedrefs[nederrors>0][i])
                    

    return nedfreqs[nederrors>0],nedfluxes[nederrors>0],nederrors[nederrors>0]
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


"""
Loading the whole catalog takes forever. If a cache exists. USE IT!
"""
if catalog.count('ned'):
    freqs,fluxes,errors = ned_spectrum(catalog,doprint=opts.v)    
else:
    specfind_cache_file = catalog[:-4]+'_'+srcname+'.vot'
    if not os.path.exists(specfind_cache_file):
        print 'specfind selection cache ',specfind_cache_file,'not found. Generating...'
        specfind = atpy.Table(catalog)
        specfind_cache = select_table_where_source(specfind,srcname)
        specfind_cache.write(specfind_cache_file)
    else:
        print "loading specfind subset cache: ",specfind_cache_file
        specfind = load_table(specfind_cache_file)    
    cat_table = select_table_where_source(specfind,srcname)
    freqs,fluxes,errors = spectrum(cat_table)

chainout = srcname+'_catalog_mcmc_chain.npz'
if not os.path.exists(chainout) or opts.mcmc_recompute:
#        print "mcmc trace file %s exists, skipping"%(chainout)
            
    print freqs,fluxes
    if len(freqs)<2:sys.exit(1)
    """
    FITTING THE MODEL
    """
    #ipdb.set_trace() 
    if model_select == 'SI': 
        print "using gaussian model"
        ndim = 3
        #set up the initial positions
        p0 = np.vstack([ #regular 2d error model
            np.random.uniform(-1,1,size=nwalkers), #alpha
            np.random.uniform(1,2,size=nwalkers),   #S0
            np.random.uniform(-1,1,size=nwalkers)
            ]).T
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI,args=[freqs,fluxes,errors])  
    
    print "burning in"
    #run the sampler 100 steps for burn-in
    pos, prob, state = sampler.run_mcmc(p0, 100)
    
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    # Starting from the final position in the burn-in chain, sample for nsteps
    # steps.
    print "sampling witn %d walkers for a %d steps"%(nwalkers,nsteps)
    sampler.run_mcmc(pos,nsteps,rstate0=state)
    print "Mean acceptance fraction:", np.mean(sampler.acceptance_fraction)
    
    alpha = sampler.flatchain[:,0]
    logS0 = sampler.flatchain[:,1]
    if n.median(logS0)>100: ipdb.set_trace()
    np.savez(chainout,catalpha=alpha,catlogS0=logS0,flatchain=sampler.flatchain)


#calculate the model parameters
trace = n.load(chainout)
logS0 = find_percentile_errors(trace['catlogS0'],confidence)
print srcname,
print a2l(n.round(10**n.array(logS0),2)),
print a2l(n.round(find_percentile_errors(trace['catalpha'],confidence),2))
print a2l(n.round(find_percentile_errors(trace['flatchain'][:,2],confidence),2))

