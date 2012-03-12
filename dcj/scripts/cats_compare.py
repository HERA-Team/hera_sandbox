#!/usr/bin/env python
#
#  cats_compare.py
#  
#
#  Created by Danny Jacobs on 12/12/09.
#  PAPER Project
#

import aipy as a, numpy as n, math as m,logging
import sys, optparse,ephem,logging,os,healpy as hpy
from pylab import * 


o = optparse.OptionParser()
a.scripting.add_standard_options(o,src=True,cal=True)
o.add_option('--sep',dest='sep',type='float',default=0.25,
    help="Maximum acceptable seperation in degrees [0.25]")
o.add_option('--helm_flags',dest='helm_flags',default=None,
    help="""Flag the catalogs based on the helmboldt meta-data with [rms cut,ncomp cut]. 
    [Default=None].  Try log10(1.25), 1.""")
o.add_option('--correct',dest='correct',default=None,
    help="Choose a catalog to be corrected against the others. [default=None]")
o.add_option('--freq_range',dest='freq_range',default="0.1,0.2",
    help="Range of frequencies (in GHz) over which to include sources. default=.100,.200")
o.add_option('--flag_outliers',dest='flag_outliers',default=None,type='float',
    help="Remove from consideration sources which differ by more than this fraction")    
o.add_option('-o',dest='outcat',
    help="Output the listed catalogs as catalog files.")
o.add_option('-p',dest='plot',action='store_true',
    help="Plot flux comparison")
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')
o.add_option('-V',dest='vverb',action='store_true',
    help='Print even more')
o.add_option('--skip_union',dest='skip_union',
    help='Skip this catalog union. eg --skip_union=nvss-vlss Can be a comma sep list.')
o.add_option('--sky',dest='sky',type='str',
    help="""Choose one catalog or comparison to plot on a mollweide projection map.
    Options are:    <catalog name> 
                    <comparison name>
                    <catalog name>_e (plots fractional error).""")
    
opts, args = o.parse_args(sys.argv[1:])

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
def union_cats(catA,catB,sep,detailed=False):
    """
    Given two catalogs of sources, return a set of paired sources that are
    seperated by sep degrees or less. Detailed also returns vector of 
    seperations.
    union = union_cats(catA,catB,sep)
    union,seps = union_cats(catA,catB,sep,detailed=True)
     """
    print len(catA),len(catB)
    union = []
    seps = []
    for k in catA.keys():
        catA_src = catA[k]
        catB_src,last_sep = find_src_in_cat(catA_src,catB,sep)
#        if k=='060633-202201': print catB_src
        if last_sep<sep * a.img.deg2rad:
            union.append([catA_src,catB_src])
            seps.append(last_sep)
    if detailed: return union,seps
    else: return union
def rms_frac_err_nx2(A):
    return n.sqrt(n.average((n.diff(A,axis=1)/A[:,0])**2))
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)
def format_src(src):
    out = src.src_name
    out += '\t' + str(src._ra) 
    out += '\t' + str(src._dec)
    out += '\t' + str(src._jys)
    return out


        

#setup logging
if opts.vverb: logging.basicConfig(level=logging.DEBUG)
elif opts.verb: logging.basicConfig(level=logging.INFO)
else: logging.basicConfig(level=logging.WARNING)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)

#add permutations to union exclusion list.
skip_union = []
if not opts.skip_union is None:
    for u in opts.skip_union.split(','):
        rev = '-'.join([u.split('-')[1],u.split('-')[0]])
        print rev
        if not rev in skip_union:
            skip_union.append(rev)
    skip_union += opts.skip_union.split(',')
else: skip_union = []



if not opts.helm_flags is None: opts.helm_flags = map(float,opts.helm_flags.split(','))
opts.freq_range = map(float,opts.freq_range.split(','))
#get the listed catalogs, if helmboldt get the spectra, etc
src_spec = a.scripting.parse_srcs(opts.src,opts.cat)
cat_names = src_spec[2]
srcs = src_spec[0]
cats = {}
#load all given catalogs
for name in cat_names:
    if not opts.cal is None:
        log.debug("looking for %s in %s.py cal file"%(name,opts.cal))
        print opts.cal,srcs,name
        cats[name]=a.cal.get_catalog(opts.cal,srcs=srcs,
            catalogs=[name])
        if len(cats[name])==0: cats[name]=a.cal.get_catalog(opts.cal,catalogs=[name])
        log.debug("found %d for %s in %s.py"%(len(cats[name]),name,opts.cal))
    else: 
        cats[name] = a.src.get_catalog(srcs=srcs,catalogs=[name])
        if len(cats[name])==0: cats[name]=a.src.get_catalog(catalogs=[name])
    log.info("loading %s witb %d entries"%( name,len(cats[name])))
for c in cats.itervalues():
    update_pos(c)
if not opts.helm_flags is None:
#load the helmboldt catalog and find the "flagged" subset
    helm_rms = opts.helm_flags[0]
    helm_ncomp = opts.helm_flags[1]
    helm = a.src.get_catalog(catalogs=['helm'])
    update_pos(helm)
    helm_meta = a._src.helm.HelmboldtCatalog()
    helm_spectrum = helm_meta.get_metadata()
    helm_rms = helm_meta.get_rms()
    helm_ncomp = helm_meta.get_ncomp()
    for name in helm.keys():
        helm[name].update_jys(0.15)
    helm_flagged = []
    for k,helm_src in helm.iteritems():
        spec = n.array(helm_spectrum[k])
        rms = helm_rms[k]
        ncomp = helm_ncomp[k]
        if rms<helm_rms and rms!=-1 and ncomp<helm_ncomp:
            helm_src.update_jys(0.15)
            helm_flagged.append(helm_src)
    helm_flagged = a._src.helm.HelmboldtCatalog(helm_flagged)   

#clean out sources not in freq range
for catname, cat in cats.iteritems():
    for name in cat.keys():
        if (cat[name].mfreq<opts.freq_range[0] 
            or cat[name].mfreq>opts.freq_range[1]):
            log.debug("%s %s %s %s %s %f"%("removing out of band source ",name,"from ", catname,\
            "at freq",cat[name].mfreq))
            null = cat.pop(name)


unions = {}
flux_unions = {}
print skip_union    
for nameA,catA in cats.iteritems():
    for nameB,catB in cats.iteritems():
        if '-'.join([nameA,nameB]) in unions.keys() or '-'.join([nameB,nameA]) in unions.keys(): continue
        #print "performing union %s U %s"%(nameA,nameB)
        if nameA!=nameB:
            union_name = '-'.join([nameA,nameB])
            if union_name in skip_union: print "skipping ",union_name; continue
            print union_name
            if not opts.helm_flags is None:
                helmf_catA_union = union_cats(helm_flagged,catA,opts.sep)
                catA = a.fit.SrcCatalog(s[1] for s in helmf_catA_union)
                update_pos(catA)
                log.info("%s %s %d"%(nameA,"(helm flagged)",len(catA)))
                helmf_catB_union = union_cats(helm_flagged,catB,opts.sep)
                catB = a.fit.SrcCatalog(s[1] for s in helmf_catB_union)
                update_pos(catB)
                print nameB,"(helm flagged)",len(catB)
            unions[union_name]=union_cats(catA,catB,opts.sep)
            if not opts.flag_outliers is None:
                temp = []
                for s in unions[union_name]:
                    if (
                    s[0]._jys/s[1]._jys>opts.flag_outliers or 
                    s[1]._jys/s[0]._jys > opts.flag_outliers): continue
                    else: temp.append([s[0]._jys,s[1]._jys])
                    print "flagging ",s[0].name,s[0]._jys, s[1].name,s[1]._jys
                flux_unions[union_name] = n.array(temp)
            else: 
                flux_unions[union_name]=n.array(
                    [[s[0]._jys,s[1]._jys] for s in unions[union_name]])

if not opts.correct is None: print "correcting ",opts.correct
print "        comparison \t m \t b \t rms error \t N"
flux_scales = {}
warnings = []
footsym = ['*','**','+','#','$']
if not opts.correct is None:
    corrected_flux_unions = {}
    for name,flux in flux_unions.iteritems():
        if name.endswith(opts.correct):
            x,y = (1,0)
            print "    ",
        elif name.startswith(opts.correct):
            x,y = (0,1)
            print "    ",
        else: 
            warn = len(warnings)
            warnings.append(footsym[warn]+" warning: finding correction for "+name.split('-')[0])
            print "%4s "%(footsym[warn],),
            x,y = (0,1)
        A = n.vstack([flux[:,x], n.ones(len(flux))]).T    
        fit = n.linalg.lstsq(A,flux[:,y])
        print "%20s: \t %2.2f \t %2.2f \t %2.2f \t\t %i" %(name,fit[0][0],fit[0][1], 
                rms_frac_err_nx2(n.vstack([flux[:,x]*fit[0][0] + fit[0][1],
                    flux[:,y]]).T), len(flux))
        flux_scales[name] = fit
        corrected_flux_unions[name] = n.vstack([flux[:,x]*fit[0][0] + fit[0][1],
                    flux[:,y]]).T
else:
    for name,flux in flux_unions.iteritems():
        print "%20s: \t %2.2f \t %2.2f \t %2.2f \t\t %i" %(name,1,0, 
                rms_frac_err_nx2(n.vstack([flux[:,0]*1 + 0,
                    flux[:,1]]).T), len(flux))
for line in warnings: print line
    


if not opts.outcat is None:
    save_cat = {}
    if opts.outcat == opts.cat:
        outfile = opts.outcat+"_cats_compare_outcat.txt"
        file = open(outfile,'w')
        for name,src in cats[opts.outcat].iteritems():
            file.write(format_src(src)+'\n')
        file.close
            
    else:
        for u in unions.keys():
            for i,cat in enumerate(u.split('-')):
                if cat==opts.outcat: save_cat[u] = i
        srcs = []
        for uname, union in unions.iteritems():
            i = save_cat[uname]
            for src_pair in union:
                srcs.append(src_pair[i])
            #cat = a.src.get_catalog(srcs=srcs,catalogs=[uname.split('-')])
            outfile = uname.split('-')[i]+'.cat.txt'
            print "writing catalog file: %s" % outfile
            file = open(outfile,'w')
            for src in srcs:
                file.write(format_src(src)+'\n')
            file.close()
if opts.plot:
    fig = figure(100)
    xlbl = []
    ylbl = []
    if not opts.correct is None:
        for name,flux in corrected_flux_unions.iteritems():
            plot(flux[:,x],flux[:,y],'.',label=name)
            xlbl.append(name.split('-')[x])
            ylbl.append(name.split('-')[y])
        xlabel(','.join(xlbl)+'[Jy]')
        ylabel(','.join(ylbl)+'[Jy]')
    else:
        for name,flux in flux_unions.iteritems():
            plot(flux[:,0],flux[:,1],'.',label=name)
            xlbl.append(name.split('-')[0])
            ylbl.append(name.split('-')[1])
        xlabel(','.join(xlbl)+'[Jy]')
        ylabel(','.join(ylbl)+'[Jy]')
    legend(loc='lower right',numpoints=1)
    #print "saving figure"
    #fig.savefig('cats_compare_1.png',dpi=600)

    draw()
    for i,name in enumerate(flux_unions):
        fig = figure(100+i+1)
        try: yerr = n.array([s[1].e_S_nu for s in unions[name]])
        except(AttributeError): yerr = None
        try: xerr =  n.array([s[0].e_S_nu for s in unions[name]])
        except(AttributeError): xerr=None
        errorbar(flux_unions[name][:,0],
            flux_unions[name][:,1],
            yerr=yerr,xerr=xerr,fmt='.')
        plot(linspace(0,n.max(flux_unions[name])),linspace(0,n.max(flux_unions[name])),'k')
        title(name)
        xlabel(name.split('-')[0])
        ylabel(name.split('-')[1])
        draw()
        if yerr is not None or xerr is not None:
            if yerr is None: 
                logL = n.average(flux_unions[name][:,0]-flux_unions[name][:,1]/ \
                    2*(xerr**2))
            elif xerr is None: 
                logL = n.average(flux_unions[name][:,0]-flux_unions[name][:,1]/ \
                    2*(yerr**2))                
            else: logL = n.average(flux_unions[name][:,0]-flux_unions[name][:,1]/ \
                2*(yerr**2 + xerr**2))
            print "logL = %07.2f"%(logL)
        del(xerr)
        del(yerr)

if not opts.sky is None:
    nside = 64
    plist = n.ones(hpy.nside2npix(nside))*0.0001
    if len(opts.sky.split('-'))==1:
        if opts.sky.endswith('_e'):
            cat = cats[opts.sky[:-2]]
            try: 
                for src in cat:
                    plist[hpy.ang2pix(nside,n.pi/2 - cat[src].dec,
                        cat[src].ra)] = cat[src].e_S_nu/cat[src]._jys
            except(AttributeError):
                log.warning("No flux errors found. Plotting catalog flux instead")
                for src in cat:
                    plist[hpy.ang2pix(nside,n.pi/2 - cat[src].dec,
                        cat[src].ra)] = cat[src]._jys
        else:
            cat = cats[opts.sky]
            for src in cat:
                plist[hpy.ang2pix(nside,n.pi/2 - cat[src].dec,
                    cat[src].ra)] = cat[src]._jys
    elif len(opts.sky.split('-'))==2:
        if opts.sky.endswith('_1'):
            for s in unions[opts.sky[:-2]]:
                src = s[1]
                ephem.FixedBody.compute(src,ephem.J2000)
                plist[hpy.ang2pix(nside,n.pi/2 - src.dec,
                    src.ra)] = src._jys
        elif opts.sky.endswith('_0'):
            for s in unions[opts.sky[:-2]]:
                src = s[0]
                ephem.FixedBody.compute(src,ephem.J2000)
                plist[hpy.ang2pix(nside,n.pi/2 - src.dec,
                    src.ra)] = src._jys
        else:
            for s in unions[opts.sky]:
                src = s[0]
                ephem.FixedBody.compute(src,ephem.J2000)
                plist[hpy.ang2pix(nside,n.pi/2 - src.dec,
                    src.ra)] = n.abs(s[0]._jys-s[1]._jys)
#    hpy.mollview(
#        n.log10(hpy.smoothing(plist,30,arcmin=1)),
#        rot=[180,45],title=opts.sky,min=10*n.min(plist))
    hpy.mollview(
        n.log10(plist),
        rot=[180,45],title=opts.sky,min=10*n.min(plist))

    show()
show()    
    
#generate downselected catalogs based on helmboldt (if requested)
#generate union catalogs for all pairs
#find flux comparison tables for each pair of catalogs
#generate/print flux-flux plots, flux corrections, rms and number of points        
    