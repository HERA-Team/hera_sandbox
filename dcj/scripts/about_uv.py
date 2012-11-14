#!/usr/bin/env python
#
#  about_uv.py
#  
#
#  Created by Danny Jacobs on 12/22/09.
#  PAPER Project
#
"""
Print a nice summary about a uv file.
"""

import aipy as a, numpy as n,math as m,os,pickle
import sys, optparse,ephem,pylab as p
from types import *
o = optparse.OptionParser()
o.set_usage("""about_uv.py [files]
Print a nice summary about a uv file""")
o.add_option('--corr_plot',action='store_true',
    help='Plot a cross-correlation matrix for the first time.')
o.add_option('--print_lst_bins',action='store_true',
    help='Print a list of the lst bins in the file [experimental]')
a.scripting.add_standard_options(o,src=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

lsts = []
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)

if not opts.src is None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.src.get_catalog(srclist,cutoff,catalogs)
    update_pos(cat)

for file in args:
    aboutfile=file+os.path.sep+'about'
    uv = a.miriad.UV(file)
    nchan = uv['nchan']
    aa = a.phs.AntennaArray((str(ephem.degrees(uv['latitud'])),
        str(ephem.degrees(uv['longitu']))),[])
    about = {'tmin':0,'tmax':0,'n_times':0,'ants':[],'pols':[]}
    dec1,dec2 = ephem.degrees(uv['dec']-30*a.img.deg2rad),\
            ephem.degrees(uv['dec']+30*a.img.deg2rad)

    if os.path.exists(aboutfile):
        about = pickle.load(open(aboutfile))
    else:
        a.scripting.uv_selector(uv,ants='auto',pol_str="yy,xx,yx,xy")
        about['tmin'] = 0
        about['tmax'] = 0
        about['n_times'] = 0
        about['ants'] = []
        about['pols'] = []
        c_time = 0

        del(uv)
        uv = a.miriad.UV(file)
        if opts.corr_plot:
            cmat = n.zeros([uv['nants']]*2)
        lsts.append([])
        for (uvw,t,(i,j)),d in uv.all():
            aa.set_jultime(t)
            lsts[-1].append(aa.sidereal_time())
            if about['tmin']==0 or t<about['tmin']: about['tmin'] = t
            if about['tmax']==0 or t>about['tmax']: about['tmax'] = t
            if c_time != t:
                c_time= t
                about['n_times'] += 1
            if not i in about['ants']:
                about['ants'].append(i)
            if a.miriad.pol2str[uv['pol']] not in about['pols']:
                about['pols'].append(a.miriad.pol2str[uv['pol']].strip())
            if opts.corr_plot and t==about['tmin']: 
                dspec = d[nchan/3:nchan*2/3]
                cmat[i,j] = n.ma.average(dspec*n.ma.conjugate(dspec))
        pickle.dump(about,open(aboutfile,'w'))
    print "-----------------------------------------"
    print file
    for k,v in about.iteritems():    
        if not type(v)==ListType: print "%10s: \t %f" % (k,v)
        else:
            print "%10s:\t"%(k,),
            for l in n.sort(v):
                print l,
            print
    print "%10s \t %3.2f MHz" % ('start freq:',uv['sfreq']*10**3,)
    print "%11s \t %3.2f MHz" %('end freq:', uv['nchan']*uv['sdf']*10**3 + uv['sfreq']*10**3,)
    print "%10s \t %3.2f kHz" % ('resolution:',uv['sdf']*10**6,)
    print "%10s \t %i " % ('nchan:',uv['nchan'])
    aa.set_jultime(about['tmin'])
    print "%10s \t %s " % ('Start LST:',ephem.hours(aa.sidereal_time()))
    print "%10s \t %s " % ('End LST:',ephem.hours(uv['lst']) )
    print "%10s \t %s \t %s" % ('dec range:',dec1,dec2)
    print "history\n",uv['history']
    print "-------------------------------------"

    #find sources in this file
    ra1,ra2 = aa.sidereal_time(),ephem.hours(uv['lst'])
    if not opts.src is None:
        print "sources in this file"
        srcs = [cat[s] for s in cat if (cat[s].ra > ra1 and cat[s].ra < ra2 \
            and cat[s].dec > dec1 and cat[s].dec < dec2)]
        if len(srcs)>0:
            for src in srcs:
                print src
    del(uv)
    if opts.corr_plot: 
        p.matshow(n.log10(cmat))
#        p.show()
if opts.print_lst_bins:
    lsts = n.diff(n.array(lsts),axis=1)
    for i in range(lsts.shape[1]):
        for j in range(lsts.shape[0]):
            print "%4.2e"%lsts[j,i],
        print 
    print lsts.shape
