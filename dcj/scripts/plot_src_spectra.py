#! /usr/bin/env python
"""
Plot the ouput of beam_src_vs_ha along with other catalog data
"""
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from capo.dcj import *
from pylab import *
from scipy import optimize
matplotlib.rcParams.update({'font.size':14})
CAT=True
ioff()
confidence=73
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
o.add_option('-v',action='store_true',
    help="Print more details about fitting process and catalog parsing")
o.add_option('--fmax',default=5e3,type='float',
    help='Maximum catalog frequency in MHz [default=5e3]')
o.add_option('--calmodel',
    help='Calibration model file.')
opts,args = o.parse_args(sys.argv[1:])
if os.path.exists(opts.calmodel):
  #get a list of the calibrators
  try:
      cals = opts.calmodel.split('gain')[0].split('_')
  except(IndexError):
      pass
def SImodel(freqs,fluxes,err=None):
    """
    input 
        frequencies, fluxes, and errors (optionaly) 
    returns 
        f_0,S0,S(f)
    if no errors are given, assumes 20% of input flux
    """
    F = n.mean(freqs)
    logS = n.log10(fluxes)
    if err is None: logerr = 0.2
    else: logerr = err/fluxes
    #fit = n.ma.polyfit(n.log10(freqs/F),logS,1)
    fitfunc = lambda p,x: p[0] + p[1] * x
    errfunc = lambda p,x,y,err: (y - fitfunc(p,x))/err
    res = optimize.leastsq(errfunc,[1.,-1.],args=(n.log10(freqs/F),logS,logerr))
    fit = res[0]
    SIM = lambda f: 10**fit[0]*(f/F)**fit[1]
    return [10**fit[0],fit[1],SIM]
def plotmodel(f0,S0,alpha):
    return lambda freq: S0*(freq/f0)**(alpha)
def find_wrapped(A,range):
    if range[0]<range[1]:
        return n.argwhere(n.logical_and(A>range[0],A<range[1]))
    else:
        return n.argwhere(n.logical_not(n.logical_and(A>range[1],A<range[0])))
def sky_sep(A,B):
    #compute distance on sphere
    #using input vectors in xyz
    (theta1,phi1) = slalib.sla_cc2s(A)
    (theta2,phi2) = slalib.sla_cc2s(B)
    return slalib.sla_sep(theta1,phi1,theta2,phi2)
def load_speclist(speclist):
    lines = open(speclist).readlines()
    speclist = {}
    for i,l in enumerate(lines):
        if l.startswith('Saving'):continue
        if l.startswith('Reading'):
            srcname = l.split('/')[-1].split('_spec')[0].strip()
            speccount=-1
            f,S,E = [],[],[]
            continue
        if l.startswith('S150'):
            S150 = float(l.split(':')[1].strip())
            continue
        if l.startswith('\\alpha'):
            SI  = float(l.split('\\alpha')[1].strip())
            speccount = i+1 
            continue
        if i==speccount:
            f.append(float(l[:3].strip()))
            S.append(float(l[9:15].strip()))
            E.append(float(l.strip()[-4:]))
            if not lines[i+1].startswith('843'):
                speccount += 1
        speclist[srcname] = {'S150':S150,'SI':SI,'freqs':f,'flux':S,'error':E}
    return speclist
def is_left(rows,cols,ind):
    A=n.zeros((rows,cols))
    A[:,0]=1
    A = A.ravel()
    return A[ind-1]
def is_bottom(rows,cols,ind):
    A=n.zeros((rows,cols))
    A[-1,:]=1
    A = A.ravel()
    return A[ind-1]
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
def filename2src(f):
    #return f.split('_')[-1]
    return f.split('_')[0].strip()
def average_spectrum(x,y,bins):
    myspectrum = []
    myerrors = []
    myfreqs = []
    for i in range(1,len(bins)):
        myspectrum.append(n.mean(y[n.logical_and(x>bins[i-1],x<bins[i])]))
        myerrors.append(n.std(y[n.logical_and(x>bins[i-1],x<bins[i])]))
        myfreqs.append((bins[i]+bins[i-1])/2)
    myspectrum=n.array(myspectrum)
    myfreqs = n.array(myfreqs)
    myerrors = n.array(myerrors)
    return myfreqs,myspectrum,myerrors
#srclist = [filename2src(f) for f in args]
#print srclist
#srclist,cutoff,catalogs, = a.scripting.parse_srcs(','.join(srclist), opts.cat)
#cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)

f=n.linspace(40,opts.fmax)
psa64freqs,PAPERcatalog = load_paper_spectra(args[0])
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
rows=5
cols=4
#load teh new points
plot_seperate = ['pic','cen']
skip = ['vela','0518-458B','1315-460','1247-401','1243-412']
for src in skip:
    try:
        PAPERcatalog.pop(src)
    except(KeyError):
        continue
multifig = 1
singlefig = 2
i=0
for srcname in sort(PAPERcatalog.keys()):
    if srcname in skip: continue
    print srcname,

    psa64fluxes,psa64errors = PAPERcatalog[srcname]['flux'],PAPERcatalog[srcname]['err']
 

    #if SED_fit.py has been run on this output already, then load the trace and plot the model possibilities
    chainin = srcname+'_mcmc_chain.npz'
    if os.path.exists(chainin):
        trace = n.load(chainin)
        logS0 = n.array(find_percentile_errors(trace['logS0'],confidence))
        alpha = n.array(find_percentile_errors(trace['alpha'],confidence))
        mcmcfreqs = n.exp(n.random.uniform(n.log(f.min()),n.log(f.max()),trace['logS0'].size))
        SED_MCMC = 10**trace['logS0']*(mcmcfreqs/150.)**trace['alpha']
        psa64_model = lambda f: 10**logS0[1]*(f/150.)**alpha[1]
    print "loading catalog data"
    if CAT:
        nedfile = srcname+'_ned_spectrum.vot'
        if os.path.exists(nedfile):
            print "plotting ned data in:",nedfile
            cat_freq,cat_flux,cat_err = ned_spectrum(nedfile,doprint=opts.v,fmax=opts.fmax*1e6)
        else:
            print "using specfind catalog"
            cat_table = select_table_where_source(specfind,srcname)
            cat_freq,cat_flux,cat_err = spectrum(cat_table)
            print len(cat_flux),len(cat_err)
    else:print
    ind = i%(rows*cols)  +1 
    fignum = int(i/(rows*cols))
    #for sources to be plotted seperately  
    if srcname in plot_seperate:
        figure(singlefig)
        clf()
        ax = subplot(111)
    else:#plot the rest in a grid plot
        figure(multifig)
        if ind==1:
            if i>0:
                figtext(0.44,0.04,'frequency [MHz]')
                figtext(0.05,0.5,'flux [Jy]',rotation='vertical')
                savefig('srcfig_%d.png'%fignum)
                savefig('srcfig_%d.eps'%fignum)
                clf()
        ax=subplot(rows,cols,ind)
    ax.set_yscale('log',nonposx='clip')
    ax.set_xscale('log',nonposy='clip')
    xlim([40,opts.fmax])
    if CAT:
        if len(cat_err)>1:
            cat_err *= 2
            cat_err_low = n.array(cat_err).copy()
            cat_err_low[(cat_flux-cat_err)<0] = cat_flux[(cat_flux-cat_err)<0]*0.99
            errorbar(cat_freq,cat_flux,yerr=[cat_err_low,cat_err],fmt='.k',capsize=0) #plot catalog data
        else:
            errorbar(cat_freq,cat_flux,yerr=cat_err,fmt='.k',capsize=0)
    psa64errors *= 2
    psa64errors_low = psa64errors.copy()
    psa64errors_low[(psa64fluxes-psa64errors)<0] = psa64fluxes[(psa64fluxes-psa64errors)<0]*0.999

    errorbar(psa64freqs,psa64fluxes,
        yerr=[psa64errors_low,psa64errors],fmt='xk',capsize=0,ms=5)
    grid()
    if os.path.exists(chainin):
        #plot(f,psa64_model(f),'k')
        plot(mcmcfreqs,SED_MCMC,',',color='0.5',alpha=0.01)
    if srcname in plot_seperate:
        xlabel('frequency [MHz]')
        ylabel('flux [Jy]')
        savefig('%s_spectrum.png'%srcname)
        #savefig('%s_spectrum.eps'%srcname)
        continue
    ylim([1,1e3])       
    #title(srcname)
    if srcname in cals:calmark='*'
    else: calmark = ''
    annotate(srcname+calmark,[0.4,0.75],xycoords='axes fraction',textcoords='axes fraction',size=12)

    #grid plot cleanup
    #for all but the bottom left, put the 10e1 (its easiest to just just bump up the bottom ylim a scosh than figure out
    #why I can't delete the goddam bottom yticklabel
    if is_left(rows,cols,ind) and not is_bottom(rows,cols,ind): ylim([1.01,1e3]) 
    tick_params(labelleft=is_left(rows,cols,ind),labelbottom=is_bottom(rows,cols,ind))

    locs,labels = yticks()
#    if not is_left(rows,cols,ind):#yticks([])    
#    if not is_bottom(rows,cols,ind) and (len(PAPERcatalog)-i)>cols:
#        xticks([])
    subplots_adjust(wspace=0,hspace=0,bottom=0.15)
    i += 1
figure(multifig)
figtext(0.44,0.04,'Frequency [MHz]')
figtext(0.05,0.5,'Flux [Jy]',rotation='vertical')
savefig('srcfig_%d.png'%(fignum+1))
#savefig('srcfig_%d.eps'%(fignum+1))





