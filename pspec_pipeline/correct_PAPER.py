#! /usr/bin/env python
"""
This is a more flexible version of the traditional correct_psa????.py script
used in compression. It does NOT overwrite metadata, unlike its predecessors
(except for antenna number & gains, which needs adjusting after removing ants).

It also assumes that there is one polarization per file (so ngains = nants)
"""

import aipy as a, numpy as n, os, sys, glob
import optparse
o = optparse.OptionParser()
#stuff you might want to use
o.add_option('-r','--rma', dest='rmants', default=[], help='Comma-separated antenna numbers to remove.')
o.add_option('-f','--flip',dest='flipants',default=[],help='Antennae that need rotating 180 degrees (data -> -1*data)')
#o.add_option('-w','--rewire',dest='rewire',default=[],help='Antennae that need pols swapped x<-->y')

#stuff you probably don't want to use
o.add_option('-j','--conj',dest='conj',action='store_true',help='Conjugate for j>i? CAUTION! This should ONLY be used on RAW data')

opts,args = o.parse_args(sys.argv[1:])

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m

rmants = map(int,opts.rmants.split(','))
if len(opts.flipants) != 0:
    flipants = map(int,opts.flipants.split(','))
else:
    flipants = []
#rewireants=map(int,opts.rewire.split(','))

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    file_jd = float('.'.join(filename.split('.')[1:-2]))
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0
    def mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        p1,p2 = a.miriad.pol2str[uv['pol']]
        
        if i in rmants or j in rmants: return p, None, None
        # prevent multiple entries arising from xy and yx on autocorrelations
        if i == j and (p1,p2) == ('y','x'): return p, None, None
        if opts.conj: #!!! CAUTION !!!
            if i > j: 
                i,j,d = j,i,n.conjugate(d)
        if i in flipants or j in flipants: d = -d  # I think the dipole for this antenna is rotated 180 deg
        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            curtime = t
        return (crd,t,(i,j)),d,f

    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'nants': uvi['nants']-len(rmants),
        'ngains': uvi['nants']-len(rmants),
        'telescop':'PAPER',
    }
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='\nCORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
