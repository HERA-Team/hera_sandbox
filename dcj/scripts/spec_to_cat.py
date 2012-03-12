#!/usr/bin/env python
#
#  spec_to_cat.py
#  
#
#  Created by Danny Jacobs on 2/17/10.
#  PAPER Project
#
"""
Collect a list of spectrum files into a catalog.
Assumes:
filename =  <prefix><srcname>.txt
non-data lines escaped by '#'
Data lines given by 
freq \t flux #units of GHz and Jys

Names correspond to entries in given catalog.

Usage:
spec_to_cat.py <files> -s all --cat=helm
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem,re
def find_in_arr(a,X):
    "Find the element of X that is closest to a. Return the index"
    try: return n.argwhere(n.abs((X-a))==n.min(n.abs(X-a))).squeeze()[0]
    except(IndexError): return n.argwhere(n.abs((X-a))==n.min(n.abs(X-a))).squeeze()

def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            del(c[s])
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

    
o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True,cal=True)
o.add_option('--prefix', dest='prefix',
    help='Filename prefix.')
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')
o.add_option('--sep',dest='sep',default=None,type='float',
    help="""If the exact source name isnt in the catalog, find the next closest within
    this radius""")
o.add_option('-m',dest='mode',default='spec',type='str',
    help="""Files can be loaded with difference parsing and file division algorithms.
    Available types: spec, channel.
    Spec: A single spectrum file for each source. [freq [GHz], flux [Jys]]
    channel: A single catalog for each frequency. [name, freq [GHz], nsamp, helm flux [Jys], pgb flux [Jys]""")
o.add_option('-f',dest='freqs',default='all',
    help="Frequency or range of frequencies or [all].")
    
opts, args = o.parse_args(sys.argv[1:])


srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)

if opts.cal != None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)
if opts.freqs=='all': freq_range = (-1,10**10)
elif opts.freqs.split('_')>1:
    freq_range = map(float,opts.freqs.split('_'))
else: freq_range = float(opts.freqs)
update_pos(cat)
print "#srcname \t ra \t dec \t flux [Jy] \t sigma [Jy] \t cat pos error [\"] \t catalog "
if opts.mode=='spec':
    for file in args:
        #open and read the spectrum
#        filelines = open(file).readlines()
#        filespec = n.array([map(float,l.split()) for l in filelines if l[0]!='#' \
#         and l[0].isdigit()]).astype(n.float)
        filespec = n.loadtxt(file)
        if len(freq_range)>1:
            chans = n.where([f>freq_range[0] and f<freq_range[1] for f in filespec[:,0]])
            filespec = filespec[chans]
        else: 
            chans = find_in_arr(freq_range,filespec[:,0])
            filespec = filespec[chans]
            filespec.shape = (1,) + filespec.shape
        #if the spectrum has error bars do a pdf computation of the median and variance
        if filespec.shape[-1]>2:
            g = lambda f,w,f0: 1/n.sqrt(2*n.pi*w**2)*n.exp(-(f-f0)**2/(2*w**2))
            flux = n.clip(n.linspace(
                n.min(filespec[:,1])-n.max(filespec[:,2]),
                n.max(filespec[:,1])+n.max(filespec[:,2]),
                num=1000),0,100000)
            pdf = n.zeros_like(flux)
            for s in filespec:
                if s[2]>0: 
                    pdf += g(flux,s[2],s[1])/n.sum(g(flux,s[2],s[1]))/n.max(n.diff(flux))
            conf_76 = n.abs(flux[n.where(
                    n.cumsum(pdf*flux**2)/n.sum(pdf*flux**2)>0.24)[0][0]] -\
                    flux[n.where(
                    n.cumsum(pdf*flux**2)/n.sum(pdf*flux**2)>(1-0.24))[0][0]])
            flux_err = conf_76
            flux = flux[find_in_arr(0.5,n.cumsum(flux*pdf)/n.sum(flux*pdf))]
            #flux_err = n.sqrt(n.median(flux**2*pdf))
        else:
            flux_err=n.sqrt(n.var(filespec[:,1]))
            flux = n.average(filespec[:,1])
        
        #compute the flux and error
        filespec[:,0] *= 1e9
        mfreq = n.average(filespec[:,0])
        #find the relevant item in the matching catalog
        if len(file.split(opts.prefix))>1:
            srcname = file.split(opts.prefix)[1].split('.txt')[0]
        else: srcname = find_srcname_in_specfile(file)
#        ra = ephem.hours(srcname.split('_')[0])
#        dec = ephem.degrees(srcname.split('_')[1])
        try: 
            src = cat[srcname]
            print "%10s \t %s \t %s \t %07.2f \t %07.2f \t 0 \t %s" % \
                (srcname,src.ra,src.dec,flux,flux_err,opts.cat)
        except(KeyError):
            if opts.verb and opts.sep is None: print srcname,"not found"
            if not opts.sep is None:
                tsrc = a.phs.RadioFixedBody(ra,dec) 
                ephem.FixedBody.compute(tsrc, ephem.J2000)
                fsrc,sep = find_src_in_cat(tsrc,cat,opts.sep)
                print "%10s \t %s \t %s \t %07.2f \t %07.2f \t %f \t %s" %\
                    (srcname,fsrc.ra,fsrc.dec,flux,flux_err,sep*180/n.pi*3600,opts.cat)
                if fsrc is None and not opts.verb is None: 
                    print "%10s not found in %s within %d deg" %\
                    (srcname,opts.cat,sep)
elif opts.mode=='channel':
    srclist={}
    for file in args:
        filelines = open(file).readlines()
        for line in filelines:
            line = line.split()
            if not srclist.has_key(line[0]): srclist[line[0]] = [(line[0],float(line[1]),float(line[4]))]
            else: srclist[line[0]].append((line[0],float(line[1]),float(line[4])))
    for src in srclist:
        s = n.array(srclist[src],dtype=[('name','S10'),('freq',n.float),('flux',n.float)])
        print "%s \t %s \t %s \t %07.2f \t %07.2f \t 0 \t %s " %\
        (s['name'][0],cat[s['name'][0]].ra,cat[s['name'][0]].dec,
            n.average(s['flux']),n.sqrt(n.var(s['flux'])),opts.cat)
        #filespec = n.array([map(float,l.split()) for l in filelines if l[0]!='#' \
#             and l[0].isdigit()]).astype(n.float)
        
                
    