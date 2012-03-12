#!/usr/bin/env python
#
#  cats_unique.py
#  
#
#  Created by Danny Jacobs on 12/12/09.
#  PAPER Project
#

import aipy as a, numpy as n, math as m
import sys, optparse,ephem
from pylab import * 
o = optparse.OptionParser()
a.scripting.add_standard_options(o,src=True)
o.add_option('--sep',dest='sep',type='float',default=0.25,
    help="Maximum acceptable seperation in degrees [0.25]")
o.add_option('--helm_flags',dest='helm_flags',default=None,
    help="""Flag the catalogs based on the helmboldt meta-data with [rms cut,ncomp cut]. 
    [Default=None].  Try log10(1.25), 1. (log10(1.25)==0.096""")
#o.add_option('--correct',dest='correct',default=None,
#    help="Choose a catalog to be corrected against the others. [default=None]")
o.add_option('--freq_range',dest='freq_range',default=".1,.2",
    help="Range of frequencies (in GHz) over which to include sources. default=.100,.200")
#o.add_option('--flag_outliers',dest='flag_outliers',default=None,type='float',
#    help="Remove from consideration sources which differ by more than this fraction")    
o.add_option('--basis',dest='basis_cat',
    help="""The catalog one wishes to use as the selection basis.  
    If not specified, the first catalog listed is used.""")
o.add_option('--list_type',dest='list_type',default='pretty',
    help="The type of source list to print. ([pretty], names)")
o.add_option('--el_cut',dest='el_cut',default="90,-90",
    help="An elevation range to search in.  eg --el_cut=90,10.  Default, 90,-90")
o.add_option('-v',dest='verb',action='store_true',
    help="print more")

opts, args = o.parse_args(sys.argv[1:])
if not opts.helm_flags is None: opts.helm_flags = map(float,opts.helm_flags.split(','))
opts.freq_range = map(float,opts.freq_range.split(','))
sep = opts.sep
opts.el_cut = map(float,opts.el_cut.split(','))
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
    if opts.verb: print len(catA),len(catB)
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
def union_compliment(catA,catB,sep):
    """
    Given two catalogs (A,B) return the subset of A for which
    there is no sources found in B.
    Input:
    catA, catB: catalog objects returned by eg a.src.get_catalog()
    sep: maximum association seperation in degrees
    Returns:
    list of sources suitable for input into a.src.get_catalog()
    """
    if opts.verb: print len(catA),len(catB)
    srcs = []
    seps = []
    for k in catA.keys():
        catA_src = catA[k]
        catB_src,last_sep = find_src_in_cat(catA_src,catB,sep)
        if catB_src is None:
            srcs.append(catA_src)
            seps.append(last_sep)
    return srcs

def rms_frac_err_nx2(A):
    return n.sqrt(n.average((n.diff(A,axis=1)/A[:,0])**2))
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)
#def keep_in_band_srcs(cat,freq_range):
#        for s in cat


#get the listed catalogs, if helmboldt get the spectra, etc
src_spec = a.scripting.parse_srcs(opts.src,opts.cat)
cat_names = src_spec[2]
srcs = src_spec[0]
cats = {}
#load all given catalogs
for name in cat_names:
    cats[name] = a.src.get_catalog(srcs=srcs,catalogs=[name])
    if len(cats[name])==0: cats[name]=a.src.get_catalog(catalogs=[name])
    if opts.verb: print name,len(cats[name])

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
            if opts.verb: print "removing out of band source ",name,"from ", catname,
            if opts.verb: print "at freq",cat[name].mfreq
            null = cat.pop(name)

        
comp_cats = {}
if not opts.basis_cat is None:
    bcat = opts.basis_cat
else:
    bcat = cat_names[0]
if opts.verb:
    if len(cats[bcat])>2e4:
        inp = raw_input("""basis catalog %s has more than 20,000 sources are 
        you sure you want to use it as a search basis? y/[n]"""%(len(cats[bcat]),))
        if not inp or inp=='n': print "wimping out!"; sys.exit()

for name,catA in cats.iteritems():
    if bcat!=name:
        srcs = union_compliment(cats[bcat],catA,sep)
        ok_srcs = [s for s in srcs if s.dec > opts.el_cut[1]*n.pi/180 and \
                    s.dec < opts.el_cut[0]*n.pi/180]

        if opts.verb: print "found %d %s sources that are not in %s, looking them up again in %s"%(
            len(ok_srcs),bcat,name,bcat)
        
        comp_cats[name] = a.src.get_catalog(catalogs=[bcat],srcs=srcs)
        if opts.list_type == "pretty":
            if opts.verb: print "sources that are in the declination range %f to %f"%(opts.el_cut[0],
                    opts.el_cut[1])
            for src in ok_srcs:
                print src
        elif opts.list_type == "names":
            for src in ok_srcs:
                print src.src_name+',',