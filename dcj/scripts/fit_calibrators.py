#! /usr/bin/env python
"""
fit a spectral model to several calibrators and a gain parameter to PAPER
used for fitting calibrator model in psa64 beamform pic strip paper

usage
catalog_spectrum_fit.py specfind_catalog.vot 2331-416

"""
import matplotlib,emcee
#matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from capo.dcj import *
from pylab import *
from scipy import optimize
import ipdb
matplotlib.rcParams.update({'font.size':18})
def idB(x):
    return 10**(x/10)
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
o.add_option('--mcmc_recompute',action='store_true',
    help='Ignore saved traces and rerun mcmc fit, otherwise just move on to plotting.')
o.add_option('--plot',action='store_true',
    help='Do MCMC contour plots')
o.add_option('-v',action='store_true',
    help="Print more details about fitting process and catalog parsing")
o.add_option('--PAPER',type='str',
    help="text file with paper spectra (output by avg_spectral_catalog.py")
o.add_option('--fmax',default=5e3,type='float',
    help="Exclude catalog values above this freq [MHz, default=5e3]")
o.add_option('--nsteps',type='int',default=1000,
    help='number of steps')
opts,args = o.parse_args(sys.argv[1:])

#Set some constants 
confidence = 73.
#set up the chain
nwalkers = 100
nsteps = opts.nsteps
nruns = 5
model_select='SI'

gridsize = 50 #number of bins in probability histogram

def lnSI_lin(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    S0 = 10**x[1]
    return  -1 * n.sum( ( fluxes - S0 * (freqs/f0)**alpha)**2 / (2 * errors**2) )
def lnSI_gain(gain,x,PAPERfreqs,PAPERfluxes,PAPERerrors):
    alpha = x[0]
    S0 = 10**x[1]
    f0 = 150.
    return -1*n.sum((PAPERfluxes*n.abs(gain) - S0 * (PAPERfreqs/f0)**alpha)**2/(2*PAPERerrors**2))
def lnSI_cal(x,catdata,PAPERdata):
    #compute the likelihood of a joint catalog fit and single gain solution for many sources
    logL = 0
    for i,srcname in enumerate(catdata):
        catprms = x[2*i:2*i+2]
        #first compute the likelihood of the catalog fit
        freqs,fluxes,errors = catdata[srcname]
        logL += lnSI_lin(catprms,freqs,fluxes,errors)
        #then compute the likelihood of the new gain parameter
        PAPERfreqs,PAPERfluxes,PAPERerrors = PAPERdata[srcname]
        logL += lnSI_gain(x[-1],catprms,PAPERfreqs,PAPERfluxes,PAPERerrors)
    return logL

    
def lnSI_exp(x,freqs,fluxes,errors):
    f0=150.
    alpha = x[0]
    logS0 = x[1]
    logS_high = n.log10(fluxes+errors)
    logS_low = n.log10((fluxes-errors).clip(.1,1e10))
    dlogS = (logS_high-logS_low)/2#the error in log space
    logS = n.log10(fluxes)
    return -1*n.sum((logS - alpha * n.log10(freqs/f0) - logS0)**2/(2*dlogS**2))#the log likelihood
lnSI = lnSI_cal
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
        spectra[srcname] = (freqs,fluxes,errors)
    return spectra

catdata = {}
for catalog in args:
    srcname = catalog.split('_')[0]
    catdata[srcname] = ned_spectrum(catalog,doprint=opts.v,fmax=opts.fmax*1e6)    
#load teh new points
PAPERdata = load_paper_spectra(opts.PAPER)
chainout = '_'.join(catdata.keys())+'_gain_mcmc_chain.npz'
if not os.path.exists(chainout) or opts.mcmc_recompute:
#        print "mcmc trace file %s exists, skipping"%(chainout)
            
    """
    FITTING THE MODEL
    """
    #ipdb.set_trace() 
    if model_select == 'SI': 
        ndim = len(catdata)*2+1
        print "fitting in %d dimensions"%(ndim)
        #set up the initial positions
        p0 = []
        for srcname in catdata:
            p0.append(np.random.uniform(-1,1,size=nwalkers)) #alpha for this src
            p0.append(np.random.uniform(1,2,size=nwalkers)) #S0 for this src
        p0.append(np.random.uniform(0,1,size=nwalkers)) #single gain parm
        p0 = np.vstack(p0).T
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnSI_cal,args=[catdata,PAPERdata])  
    
    print "burning in"
    #run the sampler 100 steps for burn-in
    pos, prob, state = sampler.run_mcmc(p0, 100)
    
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    # Starting from the final position in the burn-in chain, sample for nsteps
    # steps.
    print "sampling witn %d walkers for a %d steps"%(nwalkers,nsteps)
    for i in range(nruns):
        pos,prob,state = sampler.run_mcmc(pos,nsteps,rstate0=state)
        print "Mean acceptance fraction (run %d):"%i, np.mean(sampler.acceptance_fraction)
    
    # alpha = sampler.flatchain[:,0]
    # logS0 = sampler.flatchain[:,1]
    # if n.median(logS0)>100: ipdb.set_trace()
    np.savez(chainout,flatchain=sampler.flatchain)


#calculate the model parameters
trace = n.load(chainout)
chain = trace['flatchain']
for i,srcname in enumerate(catdata):
    print srcname,
    alpha = find_percentile_errors(chain[:,i*2],confidence)
    logS0 = n.array(find_percentile_errors(chain[:,i*2+1],confidence))
    S0 = 10**logS0
    print a2l(n.round(10**n.array(logS0),2)),
    print '%5.1f%%'%((S0[-1] - S0[0])/(2*S0[1])*100),
    print a2l(n.round(alpha,2)),
    print'%5.1f%%'%(n.abs((alpha[-1] - alpha[0])/(2*alpha[1]))*100)

print "saving chain: ",chainout
gain = n.array(find_percentile_errors(dB(n.abs(chain[:,-1])),confidence))
print "gain [dB]", a2l(n.round(gain,2))
ghist,gbins =n.histogram(dB(n.abs(chain[:,-1])),bins=100)
semilogy(gbins[1:],ghist.astype(n.float)/ghist.max(),'k')
print "histogram in ",chainout[:-4]+'_gain_conf.png'

#hist(dB(n.abs(chain[:,-1])),bins=100,histtype='step',color='k')
xlabel('flux cal [dB]')
xlim([-3,0.5])
savefig(chainout[:-4]+'_gain_conf.png')
print "gain [x]",a2l(n.round(idB(gain),2))
gain = idB(gain)
print "gain error [mult,%]", (gain[-1] - gain[0])/2/gain[1]*100

#output the newly calibrated file
PAPERcal = opts.PAPER[:-4]+'_cal.txt'
F = open(PAPERcal,'w')
print "newly calibrated catalog:",PAPERcal
F.write("#FREQS[MHz]=%s\n"%(','.join(map(str,PAPERdata[PAPERdata.keys()[0]][0]))))
for i,srcname in enumerate(sort(PAPERdata.keys())):
    freqs,fluxes,errors = PAPERdata[srcname]
    F.write("%s\t%s\t%s\n"%(srcname,
        ','.join(map(str,n.round(gain[1]*fluxes,2))),
        ','.join(map(str,n.round((gain[2]-gain[0])/2/gain[1]*fluxes+errors,3)))))
F.close()

if opts.plot:
    mcmcfreqs=n.logspace(n.log10(40),n.log10(5e6),num=chain.shape[0])
    for i,srcname in enumerate(catdata):
        figfile = srcname+'_calfit.png'
        figure(2)
        clf()
        ax = subplot(111)
        freqs,fluxes,errs = catdata[srcname]
        errorbar(freqs,fluxes,yerr=errs,fmt='.k')
        freqs,fluxes,errs = PAPERdata[srcname]
        errorbar(freqs,gain[1]*fluxes,yerr=errs,fmt='xk')
        ylim([1,500])
        xlim([40,5000])
        alpha = chain[:,i*2]
        logS0 = chain[:,i*2+1]
        SED_MCMC = 10**logS0*(mcmcfreqs/150.)**alpha
        plot(mcmcfreqs,SED_MCMC,',',color='0.5',alpha=0.01)
        ax.set_yscale('log',nonposx='clip')
        ax.set_xscale('log',nonposy='clip')
        title(srcname)
        print figfile
        savefig(figfile)
