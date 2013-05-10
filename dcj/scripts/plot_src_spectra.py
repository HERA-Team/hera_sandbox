#! /usr/bin/env python
"""
Plot the ouput of beam_src_vs_ha along with other catalog data
"""
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy,os
import optparse, sys, scipy.optimize
import capo as C
from pylab import *
from scipy import optimize

CAT=True
ioff()
confidence=73
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])
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
def select_table_where_source(table,srcname):
    #Assumes MRC names!!!
    #make a special exception for pic
    if srcname=='pic':
        seq=20658
    elif srcname=='cen':
        seq=61767
    elif srcname=='hyd':
        seq=43497
    else:
        result = table.where(table['Name']== 'MRC %s'%srcname)
        seq = result['Seq']
    return table.where(table['Seq']==seq)
def spectrum(table):
    return table['nu'],table['S_nu_']/1e3,table['e_S_nu_']/1e3
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
#srclist = [filename2src(f) for f in args]
#print srclist
#srclist,cutoff,catalogs, = a.scripting.parse_srcs(','.join(srclist), opts.cat)
#cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)

f=n.linspace(40,2e3)
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
plot_seperate = ['0518-458B','pic','cen']
skip = ['vela']
#try:
#    PAPERcatalog.pop('0518-458B')
#    PAPERcatalog.pop('pic')
#    PAPERcatalog.pop('cen')
#    PAPERcatalog.pop('vela')
#except(KeyError):
#    pass
multifig = 1
singlefig = 2
i=0
for srcname in sort(PAPERcatalog.keys()):
    if srcname in skip: continue
    print srcname,
#    srcname=filename2src(filename)
#    D = n.load(filename)
#    freq = D['freq']*1e3
#    spec = n.ma.masked_invalid(D['spec'])
#    psa64freqs,psa64fluxes,psa64errors = average_spectrum(freq,spec,freqbins)
    psa64fluxes,psa64errors = PAPERcatalog[srcname]['flux'],PAPERcatalog[srcname]['err']
 
    
#    if i==0:MYCATALOG.write("#FREQS[MHz]=%s\n"%(','.join(map(str,psa64freqs))))
#    #srcflux,srcSI,specmodel = SImodel(freq,n.real(spec),n.imag(spec)) #fit a spectral model
#    #psa64_S,psa64_SI,psa64_model = SImodel(psa64freqs,n.real(psa64fluxes),psa64errors)
#    MYCATALOG.write("%s\t%s\t%s\n"%(srcname,
#        ','.join(map(str,n.round(n.real(psa64fluxes),2))),
#        ','.join(map(str,n.round(psa64errors,3)))))
#
    #if SED_fit.py has been run on this output already, then load the trace and plot the model possibilities
    chainin = srcname+'_mcmc_chain.npz'
    if os.path.exists(chainin):
        trace = n.load(chainin)
        logS0 = n.array(find_percentile_errors(trace['logS0'],confidence))
        alpha = n.array(find_percentile_errors(trace['alpha'],confidence))
        mcmcfreqs = n.exp(n.random.uniform(n.log(f.min()),n.log(f.max()),trace['logS0'].size))
        SED_MCMC = 10**trace['logS0']*(mcmcfreqs/150.)**trace['alpha']
        psa64_model = lambda f: 10**logS0[1]*(f/150.)**alpha[1]

    if CAT:
        if srcname!='pic':
            print "using specfind catalog"
            cat_table = select_table_where_source(specfind,srcname)
            cat_freq,cat_flux,cat_err = spectrum(cat_table)
        else:
            print "using NED data"
            cat_freq,cat_flux,cat_err = pic_spectrum('pic_ned_spectrum.xml')    
        try: 
            cat_flux[0]
            catS,catalpha,catmodel    = SImodel(cat_freq,cat_flux,cat_err)
        except(IndexError): catS,catalpha,catmodel = (cat_flux,-0.85,lambda x: cat_flux*(x/cat_freq/1e3)**-0.85)#by default use a power law
        print n.round(catmodel(150.),2),n.round(catalpha,3)
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
    xlim([40,2000])
     
    if CAT: errorbar(cat_freq,cat_flux,yerr=cat_err,fmt='.') #plot catalog data
    psa64errors_low = psa64errors.copy()
    psa64errors_low[(psa64fluxes-psa64errors)<0] = psa64fluxes[(psa64fluxes-psa64errors)<0]*0.999
    errorbar(psa64freqs,psa64fluxes,
        yerr=[psa64errors_low,psa64errors],fmt='.')
    if os.path.exists(chainin):
        plot(f,psa64_model(f),'k')
        plot(mcmcfreqs,SED_MCMC,',',color='b',alpha=0.01)
    if srcname in plot_seperate:
        xlabel('frequency [MHz]')
        ylabel('flux [Jy]')
        savefig('%s_spectrum.png'%srcname)
        #savefig('%s_spectrum.eps'%srcname)
        continue
    ylim([1,1e3])

    annotate(srcname,[0.5,0.75],xycoords='axes fraction',textcoords='axes fraction',size=9)
    #grid plot cleanup
    #for all but the bottom left, put the 10e1 (its easiest to just just bump up the bottom ylim a scosh than figure out
    #why I can't delete the goddam bottom yticklabel
    if is_left(rows,cols,ind) and not is_bottom(rows,cols,ind): ylim([1.01,1e3]) 
    locs,labels = yticks()
    if not is_left(rows,cols,ind):yticks([])
    if not is_bottom(rows,cols,ind) and (len(PAPERcatalog)-i)>cols:
        xticks([])
    subplots_adjust(wspace=0,hspace=0)
    i += 1
figure(multifig)
figtext(0.44,0.04,'frequency [MHz]')
figtext(0.05,0.5,'flux [Jy]',rotation='vertical')
savefig('srcfig_%d.png'%(fignum+1))
#savefig('srcfig_%d.eps'%(fignum+1))





