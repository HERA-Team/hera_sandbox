#! /usr/bin/env python
import aipy as a, numpy as np
import sys,optparse,ephem
from capo import arp,dfm

#DEBUG = True
DEBUG = False
if DEBUG: from pylab import *

o = optparse.OptionParser()
o.set_usage('add_lst_1x.py -C [cal] --lst_res=[s] files')
a.scripting.add_standard_options(o,cal=True)
o.add_option('--lst_res',dest='dlst',type='float',default=30.,help='Resolution in seconds of output LST bins (default 30)')
o.add_option('--lst_rng', default='0_24', help='Range of LSTs to bin (in hours).')
o.add_option('--tfile',dest='tfile',type='float',default=600.,help='Length of time spanned by each input file.')
opts,args = o.parse_args(sys.argv[1:])

#Create AntennaArray
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
del(uv)

ra1,ra2 = map(ephem.hours,opts.lst_rng.split('_'))

def in_lst_range(lst):
    if ra1 < ra2: return lst >= ra1 and lst < ra2
    else: return lst >= ra1 or lst < ra2

print 'Filtering input files for LSTs of interest'
nargs = []
for f in args:
    uv = a.miriad.UV(f)
    (crd,t,bl),_d,_f = uv.read(raw=True)
    aa.set_jultime(t)
    start_t = aa.sidereal_time()
    aa.set_jultime(t + opts.tfile*a.ephem.second)
    end_t = aa.sidereal_time()
    if start_t < end_t:
        if ra1 < ra2:
            if end_t < ra1 or start_t > ra2: continue
        else:
            if end_t < ra1 and start_t > ra2: continue
    else:
        if ra1 < ra2:
            if start_t > ra2 and end_t < ra1: continue
    nargs.append(f)

#Read in Data
crds,Sx1,Sx2,wgt = {},{},{},{}
for file in nargs:
    uv = a.miriad.UV(file)
    if DEBUG: a.scripting.uv_selector(uv,'0_8','xx')
    print 'Reading',file
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        aa.set_jultime(t)
        jd = int(t)
        lst = aa.sidereal_time()
        bin = dfm.lst2bin(lst,bin_width=opts.dlst)
        if not in_lst_range(lst): continue
        blp = a.pol.ijp2blp(i,j,uv['pol'])
        if not blp in Sx1.keys():
            Sx1[blp] = {}
            Sx2[blp] = {}
            wgt[blp] = {}
        if not bin in Sx1[blp].keys():
            Sx1[blp][bin] = 0. 
            Sx2[blp][bin] = 0. 
            wgt[blp][bin] = 0. 
        d = np.where(d,d,0.)
        Sx1[blp][bin] += d 
        Sx2[blp][bin] += np.abs(d)**2
        wgt[blp][bin] += dfm.bit_flip(f) 
        crds[blp] = uvw

outfiles = ['sum','sumsq','wgt']

uvi = a.miriad.UV(nargs[0])
for of in outfiles:
    filename = 'lst.%s.%s.uvT'%(of,opts.lst_rng)
    print '\t Writing',filename
    data = {'sum':Sx1,'sumsq':Sx2,'wgt':wgt}[of]
    uvo = a.miriad.UV(filename,status='new')
    uvo.init_from_uv(uvi)
    for blp in data:
        i,j,uvo['pol'] = a.pol.blp2ijp(blp)
        lstbins = data[blp].keys()
        lstbins.sort()
        for l in lstbins:
            uvo['lst'],uvo['ra'],uvo['obsra'] = lst,lst,lst
            preamble = (crds[blp],l,(i,j))
            uvo.write(preamble,data[blp][l],np.zeros(data[blp][l].shape))

if DEBUG:
    chan = 85
    blp = Sx1.keys()[0]
    figure()
    #subplot(211)
    s1 = np.array([Sx1[blp][lst] for lst in lstbins])
    s2 = np.array([Sx2[blp][lst] for lst in lstbins])
    _w = np.array([wgt[blp][lst] for lst in lstbins])
    mean = s1/_w
    var  = s2/_w
    var -= np.abs(mean)**2

    for jd in LST.keys(): 
        x = LST[jd]
        x.sort()
        y = np.array([np.abs(D[jd][blp][lst][chan]) for lst in x])
        plot(x,y,'.')
    
    plot(lstbins,np.abs(mean[:,chan]),'k-')
    plot(lstbins,np.abs(mean[:,chan]) - np.sqrt(var)[:,chan],'k:')
    plot(lstbins,np.abs(mean[:,chan]) + np.sqrt(var)[:,chan],'k:')
   
    ylabel('|V| (Jy)')

    xlabel('LST bin number (binned every %d seconds, starting at noon' % opts.dlst)

    draw()

if DEBUG: 
    #savefig('training_set.eps',fmt='eps')
    show()
