#python script for plotting the pictor A spectrum figure in psa64 spectra pic stripe paper


import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from capo.dcj import *
from pylab import *
from scipy import optimize
matplotlib.rcParams.update({'font.size':14})
def load_paper_spectra(filename):
    lines = open(filename).readlines()
    spectra ={}
    for i,line in enumerate(lines):
        try:
            if line.startswith('#'):
                freqs = n.array(map(float,line.split('=')[1].split(','))).squeeze()
                continue
            line = line.split()
            srcname = line[0]
            fluxes = n.array(map(float,line[1].split(',')))
            errors = n.array(map(float,line[2].split(',')))
            spectra[srcname] = {'flux':fluxes,'err':errors}
        except:
            print "error on line ",i,':',line
            raise
    return freqs,spectra
o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('-v',action='store_true',
    help="Print more details about fitting process and catalog parsing")
o.add_option('--fmax',default=5e3,type='float',
    help='Maximum catalog frequency in MHz [default=5e3]')
o.add_option('--gain',default=1.0,type='float',
    help='gain to apply to raw pictor spectrum')
opts,args = o.parse_args(sys.argv[1:])    

psa64freqs,PAPERcatalog = load_paper_spectra(args[0])
srcname='pic'
psa64fluxes,psa64errors = PAPERcatalog[srcname]['flux'],PAPERcatalog[srcname]['err']
confidence = 76
f=n.linspace(40,opts.fmax)
#figure(figsize=(16,12))
ax = subplot(111)
chainin = srcname+'_mcmc_chain.npz'
if os.path.exists(chainin):
    trace = n.load(chainin)
    logS0 = n.array(find_percentile_errors(trace['logS0'],confidence))
    alpha = n.array(find_percentile_errors(trace['alpha'],confidence))
    mcmcfreqs = n.exp(n.random.uniform(n.log(f.min()),n.log(f.max()),trace['logS0'].size))
    SED_MCMC = 10**trace['logS0']*(mcmcfreqs/150.)**trace['alpha']
    psa64_model = lambda f: 10**logS0[1]*(f/150.)**alpha[1]
    plot(mcmcfreqs,SED_MCMC,',',color='0.5',alpha=0.01)
nedfile = srcname+'_ned_spectrum.vot'
print "plotting ned data in:",nedfile
cat_freq,cat_flux,cat_err = ned_spectrum(nedfile,doprint=opts.v,fmax=opts.fmax*1e6)

cat_err *= 2
cat_err_low = n.array(cat_err).copy()
cat_err_low[(cat_flux-cat_err)<0] = cat_flux[(cat_flux-cat_err)<0]*0.99
errorbar(cat_freq,cat_flux,yerr=[cat_err_low,cat_err],fmt='.k',capsize=0) #plot catalog dat

psa64errors *= 2
psa64errors_low = psa64errors.copy()
psa64errors_low[(psa64fluxes-psa64errors)<0] = psa64fluxes[(psa64fluxes-psa64errors)<0]*0.999
errorbar(psa64freqs,psa64fluxes,
    yerr=[psa64errors_low,psa64errors],fmt='xk',capsize=0,ms=5)
xlabel('frequency [MHz]')
ylabel('flux [Jy]')
ax.set_yscale('log',nonposx='clip')
ax.set_xscale('log',nonposy='clip')
ylim([9,1e3])
xlim([40,opts.fmax])
savefig('pictor_spectrum.png')
D = n.load('pic_spec.npz')
freqs,spec = D['freq'],D['spec']
print freqs[0]
if True:
    clf()
    inset = subplot(111)
else:
    inset=axes([0.18,0.2,0.3,0.3])

#inset.plot(freqs*1e3,spec)
chans = n.argwhere(n.logical_and(freqs>0.125,freqs<0.170))
inset.errorbar(cat_freq,cat_flux,yerr=[cat_err_low,cat_err],fmt='.k',capsize=0)
inset.semilogy(freqs[chans]*1e3,spec[chans]*opts.gain,'.',color='0.5')
inset.errorbar(psa64freqs,psa64fluxes,
    yerr=[psa64errors_low,psa64errors],fmt='xk',capsize=0,ms=8)
inset.plot(freqs*1e3,psa64_model(freqs*1e3),'k')
inset.set_ylim([90,1e3])
inset.set_xlim([100,200])
#inset.set_yticks([300,400,500,600])
inset.set_xticks([110,130,150,170,190])
print inset.get_yticklabels()
#inset.set_yticklabels(map(str,[300,400,500,600]))
xlabel('Freq [MHz]')
ylabel('Flux [Jy]')
savefig('pictor_spectrum_zoom.png')


