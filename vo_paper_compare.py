#!/usr/bin/env python
#
#  vo_paper_compare.py
#  
#
#  Created by Danny Jacobs on 2/18/10.
#  PAPER Project
#
"""
#plot catalog values along with any data found in spectra files or the given
#vo URI

VO desktop must be running to access VO data.

NB: spectra files not yet supported
"""
import aipy as a, numpy as n,math as m,time,logging,os,atpy
import sys, optparse,vo.table as vot,ephem as ep,re,warnings
from pylab import *
from astrogrid import ConeSearch
from cStringIO import StringIO
from numpy import *
from numpy import log as loge
warnings.simplefilter('ignore',n.RankWarning)
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('--sep',dest='sep',type='float',default=0.1,
    help="Maximum acceptable seperation in degrees [0.1]")
o.add_option('--uri',default="ivo://CDS.VizieR/VIII/85",
    help="URI for vo data source that supports cone search. default=ivo://CDS.VizieR/VIII/85 [specfindv2]")
o.add_option('--fig',dest='fig',type='int',default=32,
    help="Figure number at which to start [32]")
o.add_option('--max_sub',dest='max_sub',type='int',default=25,
    help="Max number of subplots per figure [25]")
o.add_option('--prefix',dest='prefix',type='str',
    help="Spectrum file prefix default=<cal file name>. If that doesn't work, looks for a name in the file.")
o.add_option('-p',dest='plot',action='store_true',
    help="Plot.")
o.add_option('--fmax',type='float',
    help="Maximum frequency to search for in the vo. [GHz]")
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")
o.add_option('--res',type='float',
    help="Sky resolution in degrees. Should be close to expected source density. [auto]")
o.add_option('--time',type='float',
    help="Time at which to calculate catalog. [jd]  Adds 'sun_sep' and 'transit' to vot output.")
opts, args = o.parse_args(sys.argv[1:])


def unit(votarray,item):
    units = {'mJys':1/1000.,
            'mJy':1/1000.,
            'MHz':1/1000.}
    for i,t in enumerate(votarray.array.dtype.names):
        if t.startswith(item):uname=votarray.fields[i].unit
    try: return units[uname]
    except(KeyError): 1.
def find_src_name(names):
    name_prefs = ('3C','VLSS','NVSS','TXS')
    for cat_name in name_prefs:
        for name in names:
            if name.startswith(cat_name): return name
    return None
def find_src_in_cat(src,cat,sep):
    """
    Given a src (src = aipy.amp.RadioFixedObject or lower) 
    and a catalog (cat = aipy.amp.RadioCatalog or lower), find the closest
    source in cat that falls within sep [deg].
    Return: closest source object and seperation (ephem.Angle).
    If no source found within sep, return None,closest_sep
    """
    last_sep=n.pi
    closest_src=None
    for nk,test_src in cat.iteritems():
        cursep = ephem.separation(src, test_src)
        if cursep <= last_sep:
            last_sep = cursep
            closest_src = test_src
    if last_sep<sep * a.img.deg2rad: return closest_src,last_sep
    else: return None,last_sep
def find_srcname_in_specfile(filename):
    for l in open(filename).readlines():
        m = re.search("name=([a-zA-Z0-9_+]+)",l)
        if not m is None: 
            return m.groups()[0]
def find_seq_of_brightest_src(vodata):
    return vodata[n.where(vodata['S_nu_']==n.max(vodata['S_nu_']))]['Seq'][0]

if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)

    
if opts.prefix is None and not opts.cal is None:
    opts.prefix=opts.cal+'_'
elif opts.prefix is None:
    opts.prefix = 'none'
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
#print srclist,cutoff,catalogs
if opts.cal != None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)
#print dir(cat[cat.keys()[0]])
#print [cat[src]._jys for src in cat]
#try: print [cat[src].e_S_nu for src in cat]
#except(AttributeError): print "this catalog does not have flux errors"
#load any spectra files and associated spectra with catalog sources

file_srcs = []
spectra = {}
for file in args:
    if len(file.split(opts.prefix))>1:
        srcname = file.split(opts.prefix)[1].split('.txt')[0]
    else: srcname = find_srcname_in_specfile(file)
    if srcname.startswith('3c'):srcname = srcname[2:]
#        ra = srcname.split('_')[0]
#        dec = srcname.split('_')[1]
#        src = a.phs.RadioFixedBody(ra,dec,srcname=srcname)
#        src,sep = find_src_in_cat(src,cat,opts.sep)
#        file_srcs.append(src)
    spectra[srcname] = n.loadtxt(file)
    #n.array(
#        [(float(s.split()[0]),float(s.split()[1])) for s in open(file).readlines() if s[0]!='#' \
#         and s[0].isdigit()]).astype(n.float)
#print spectra
#    scat = a.phs.RadioCatalog(file_srcs)
curfig = opts.fig
if opts.plot: 
    fig = figure(curfig,figsize=(15,15))
    clf()
nsub = n.min((opts.max_sub,len(cat)))
m1 = round(sqrt(nsub))
m2 = ceil(nsub/m1)
#print m1,m2
print '\t'.join(('Name','Ra','Ra_rad','Dec','Dec_rad',
                'S_sf','e_S_sf','S_nu_paper','e_S_nu_paper',
                'logL_spec','logL_cat','seq','search_r','sep',
                'altcat','altname','n_Seq',
                'col','row','fig'))
Nun = "\'\'"
newcols = [
            {'name':'PAPER_Name','units':'','dtype':'S20'},
            {'name':'PAPER_Ra','units':'radians','dtype':float},
            {'name':'PAPER_Dec','units':'radians','dtype':float},
            {'name':'S_sf','units':'Jys','dtype':float},
            {'name':'e_S_sf','units':'Jys','dtype':float},
            {'name':'S_nu_paper','units':'Jys','dtype':float},
            {'name':'e_S_nu_paper','units':'Jys','dtype':float},
            {'name':'logL_spec','units':'likely','dtype':float},
            {'name':'logL_cat','units':'likely','dtype':float},
            {'name':'search_r','units':'degrees','dtype':float},
            {'name':'col','units':'','dtype':int},
            {'name':'row','units':'degrees','dtype':int},            
            {'name':'fig','units':'degrees','dtype':int},
            {'name':'PAPER_seq','units':'','dtype':int},
            {'name':'sf_select','units':'logic','dtype':bool}
            ]
if not opts.time is None:
    newcols.append({'name':'transit','units':'hour','dtype':float})
    newcols.append({'name':'sun_sep','units':'deg','dtype':float})
    sun =ep.Sun()
    aa = a.cal.get_aa(opts.cal,0.1,0.1,1)
    aa.set_jultime(opts.time)
    sun.compute(aa)
try: 
    notfound = []
    for i,src in enumerate(cat):
        catstr = ''
        src=cat[src]
        ep.FixedBody.compute(src,ep.J2000)
        src.mfreq = 0.15
        src.update_jys(0.15)
        if mod(i,opts.max_sub)==0 and i!=0: 
            curfig +=1
            if opts.plot:
                figtext(0.05,0.5,"flux [Jys]",rotation='vertical')
                figtext(0.5,0.05,"frequency [MHz]")
                subplots_adjust(hspace=0,wspace=0)
                draw()
                fig.savefig(opts.cat+"_vo_paper_compare_"+str(curfig)+".png",dpi=500)
                fig=figure(curfig,figsize=(15,15))
                clf()
        catstr += src.src_name+"\t"
        if opts.plot:
            ax = subplot(m2,m1,i%nsub+1,xscale='log',yscale='log')
            try: errorbar(src.mfreq,src.jys,yerr=src.e_S_nu,fmt='.',linewidth=3)
            except(AttributeError): 
             #   print "This catalog does not have flux errors"
                plot(src.mfreq,src._jys)
            if spectra.has_key(src.src_name):
                spectrum = spectra[src.src_name]
                if spectrum.shape[1]==3:
                    errorbar(spectrum[:,0]*10**3,spectrum[:,1],yerr=spectrum[:,2],
                    fmt='.',alpha=0.3)
                else:
                    plot(spectrum[:,0]*10**3,spectrum[:,1],'.')
        if not opts.uri is None:
         #   print "performing cone search"
            catstr += "%s \t %8.5f \t %s \t %8.5f \t"%(src.ra, float(src.ra),src.dec, float(src.dec))
            #perform the cone search and make a spectrum numpy array
            log.debug("Connecting to %s"%(opts.uri))
            conesrch = ConeSearch(opts.uri)
            log.debug("Performing conesearch")
            votable = conesrch.execute(src.ra/ep.degree,src.dec/ep.degree,opts.sep)
            log.debug("Converting votable string to file object")
            votable_f = StringIO(votable)
            log.debug("Parsing table from votable to numpy.")
            try: T
            except(NameError):
                T = atpy.Table()
                T.vo_read(votable_f,pedantic=False)
                votable_f.reset()
            VO = atpy.Table()
            VO.vo_read(votable_f,pedantic=False)
            if len(VO)==0: notfound.append(catstr); continue
            for col in newcols:
                VO.add_empty_column(col['name'],col['dtype'],unit=col['units'],null=Nun)
                if not col['name'] in T.data.dtype.names:
                    T.add_empty_column(col['name'],col['dtype'],unit=col['units'],null=Nun)
            votable_f.reset()
            votable_np = vot.parse_single_table(votable_f,pedantic=False)
            
            vodata = votable_np.array
            Jys = unit(votable_np,'S_nu_')
            eJys = unit(votable_np,'e_S_nu_')
            MHz = unit(votable_np,'nu')       
            if len(vodata)<2: notfound.append(catstr); continue
            #extract a single spectrum (hope there was only one hit, throw a warning at the end if not)
            seqs = list(set(vodata['Seq']))
            vodata = sort(vodata[n.where(vodata['nu']*MHz<opts.fmax)[0]],order='nu')
            #compute best fits for each selected source
            P = {}
            res = {}
            ill = []
            e_P = {}
            for seq in seqs: #calculate the fit probability for each source found
                spec  = vodata[n.where(vodata['Seq']==seq)[0]]
                if len(spec)>2:
#                    print spec['nu']*MHz,spec['S_nu_']*Jys,
                    out = n.polyfit(n.log10(spec['nu']*MHz),n.log10(spec['S_nu_']*Jys),1,full=True)
                    P[seq] = n.poly1d(out[0])
                    res[seq] = sqrt(out[1]/len(spec['nu']))
                    e_P[seq] = lambda f: n.max(10**(P[seq](log10(f))+res[seq]) -\
                                10**P[seq](log10(f)))
                else: 
                    ill.append(seq)
            if len(ill)==len(seqs):
                print '#',catstr
            #compute likelihoods for sources
            seqs = [seq for seq in seqs if not seq in ill]
            if opts.res is None: sres = 1./sqrt(len(seqs)/(pi*opts.sep**2))

            try: 
                src_logLs = []
                for seq in seqs:
                    spec = n.sort(vodata[n.where(vodata['Seq']==seq)[0]],order='nu').squeeze()
                    src_logL = -(log10(src.jys)-P[seq](log10(src.mfreq)))**2/\
                         (2*((log10(src.e_S_nu+src.jys)-log10(src.jys))**2+e_P[seq](src.mfreq)**2)) +\
                         -n.max(spec['_r'])**2/(2*sres**2)
                    src_logLs.append(src_logL)
                src_logLs =n.array(src_logLs)
                src_logLs -= loge(sum(exp(src_logLs)))
                for seq in seqs:
                    VO.data['logL_cat'][n.where(VO.data['Seq']==seq)[0]] = src_logLs
                    VO.data['S_sf'][n.where(VO.data['Seq']==seq)[0]]=10**P[seq](log10(src.mfreq))
                    VO.data['e_S_sf'][n.where(VO.data['Seq']==seq)[0]]=\
                            10**(P[seq](log10(src.mfreq))+res[seq]) -\
                                10**P[seq](log10(src.mfreq))
                    VO.data['S_nu_paper'][n.where(VO.data['Seq']==seq)[0]] = src.jys
                    VO.data['e_S_nu_paper'][n.where(VO.data['Seq']==seq)[0]] = src.e_S_nu
                src_logL = n.max(src_logLs)                    
                seq = seqs[n.where(src_logLs==src_logL)[0]]
                VO.data['PAPER_seq'] = seq
                VO.data['PAPER_Ra'][n.where(VO.data['Seq']==seq)[0]] = repr(src.ra)
                VO.data['PAPER_Dec'][n.where(VO.data['Seq']==seq)[0]] = repr(src.dec)
                VO.data['PAPER_Name'][n.where(VO.data['Seq']==seq)[0]] = src.src_name
                VO.data['sf_select'][n.where(VO.data['Seq']==seq)[0]] = 1
                VO.data['search_r'] = sres
#                if seq==10068: 
#                    print vodata[n.where(vodata['Seq']==seq)[0]]
#                    print '.'*10
#                    print VO.data
            except(AttributeError): 
             #   print "This catalog does not have flux errors"
                src_logL = -99
                #if not enough flux info, use the nearest position only
                seq = vodata[n.where(vodata['_r']==n.min(vodata['_r']))[0]]['Seq'][0]    
            catstr += "%6.3f\t"%(10**P[seq](log10(src.mfreq)))
            catstr += "%6.4f\t"%(10**(P[seq](log10(src.mfreq))+res[seq]) -\
                10**P[seq](log10(src.mfreq)))# n.abs(10**P[seq](log10(src.mfreq)) - \
#                    10**P[seq](log10(src.mfreq)-res[seq])))#-res[seq])
            catstr += "%6.3f\t"%(src.jys)
            catstr += "%6.4f\t"%(src.e_S_nu)            
            if spectra.has_key(src.src_name):
                spectrum = spectra[src.src_name]
                spectrum = spectrum[n.where(spectrum[:,2]>0)[0]] 
                if spectrum.shape[-1]>2:
                    spec_logL = n.sum(-(n.log10(spectrum[:,0])- \
                        P[seq](n.log10(spectrum[:,0])))**2/(2*log10((spectrum[:0,]-spectrum[:,2]))**2))
                    spec_logL += -loge(exp(-spec_logL)+sum(Prob.values()))
                    VO.data['logL_spec'][n.where(vodata['Seq']==seq)[0]] = spec_logL
                    catstr += spec_logL+"\t"
                    #print n.average(spectrum[:,2]/spectrum[:,1]),"\t",
            else: catstr += "\'\'\t"       
            if n.isinf(src_logL):
                catstr += "0"+"\t"
            else:
                catstr += str(src_logL)+"\t"
            catstr += str(seq)+"\t"
            catstr += str(sres)+ "\t"
            catstr += str(max(spec['_r'])) + "\t"
            #the rest takes care of the nice subplot with labels
            if len(vodata['S_nu_'])<1 and opts.plot: 
                print "no data!", 
                ylim([0.1,500])
                xlim([0.050,opts.fmax+1000])
                text(100,0.2,src.src_name,size='small',weight='normal')
                if (i % m1): yticks([])
                if ((m2-1)*m1>(i%opts.max_sub+(m2*m1-nsub))): xticks([])
                continue
            spec = n.sort(vodata[n.argwhere(vodata['Seq']==seq)],order='nu').squeeze()
            if vodata['n_Seq'][0]==0: fmt = '.'; ms =7
            else: fmt='s';ms=5
            if opts.plot:
                errorbar(spec['nu']*MHz,spec['S_nu_']*Jys,yerr=spec['e_S_nu_']*eJys\
                        ,fmt=fmt,markersize=ms,linewidth=2)
                plot(linspace(0.050,opts.fmax+1000),
                    10**(P[seq](log10(linspace(0.050,opts.fmax+1000)))),'-k')
                ylim([0.1,10000])
                xlim([0.050,opts.fmax+1000])
            #apply src label                
            srcname = find_src_name(spec['Name'])
            if srcname is None: srcname = 'PAPER '+ src.src_name
            catstr += srcname+'\t'+str(n.max(spec['n_Seq']))+'\t'
            catstr += str(i % m1)+'\t'+str(int(i)/int(m1)) + '\t' + str(curfig)
            VO.data['col'][n.where(VO.data['Seq']==seq)[0]] = i%m1+1
            VO.data['row'][n.where(VO.data['Seq']==seq)[0]] = int(i)/int(m1)+1
            VO.data['fig'][n.where(VO.data['Seq']==seq)[0]] =curfig
            if not opts.time is None:
                src.compute(aa)
                VO.data['sun_sep'] = ep.separation(src,sun)
                VO.data['transit'] = aa.next_transit(src).tuple()[3]
            if opts.plot:
                text(0.1,0.85,srcname,size=8,weight='normal',transform=ax.transAxes)
        if (i % m1) and opts.plot: yticks([])
        if ((m2-1)*m1>(i%nsub+(m2*m1-nsub))) and opts.plot: xticks([])
        if opts.plot: draw()
        print catstr
#        print VO.data['Seq'],VO.data['logL_cat']
        T.append(VO)
except(KeyboardInterrupt): sys.exit()

T.write('test.vot')
print '# sources with no match'
for catstr in notfound:
    print '#'+catstr
#subplots_adjust(hspace=0.35,wspace=0.25)
if opts.plot:
    figtext(0.05,0.5,"flux [Jys]",rotation='vertical')
    figtext(0.5,0.05,"frequency [MHz]")
    subplots_adjust(hspace=0,wspace=0)
    draw()
    fig.savefig(opts.cat+"_vo_paper_compare_"+str(curfig)+".png",dpi=500)
    show()