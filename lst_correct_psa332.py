#!/usr/bin/env python
#
#  lst_correct_psa332.py
#  
#
#  Created by Danny Jacobs on 5/20/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
o.set_usage('lst_correct_psa332.py <files>')
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])
aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9'))
curtime = 0
for filename in args:
    print filename, '->', filename+'l'
    if os.path.exists(filename+'l'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'l', status='new') 
    
    def mfunc(uv,p,d):
        crd,t,(i,j) = p
        if t!= curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
        return p,d
    override={'telescope':'PAPER'}
    uvo.init_from_uv(uvi,override=override)
    uvo.pipe(uvi,mfunc=mfunc,
        append2hist='LST_CORRECT_PSA332:   Fixed LST/RA/obsra.')
    del(uvo)