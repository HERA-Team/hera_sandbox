#!/usr/bin/env python
#
#  second_correct_psa331.py
#  
#
#  Created by Danny Jacobs on 6/23/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'C'
    if os.path.exists(filename+'C'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
#    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
#    flagants = map(a.miriad.bl2ij,[t[0] for t in a.scripting.parse_ants(opts.ant,uvi['nants'])])
#    chan = a.scripting.parse_chans(opts.chan,uvi['nchan'])
    uvo = a.miriad.UV(filename+'C', status='new')
    def mfunc(uv, p, d):
        uvw,t,(i,j) = p
        if (j-i)>=18: 
            d = n.conjugate(d); 
            print i,j,i-j,'*'
            if i%2 != j%2:
                print i,j,'-->',                    
                i += (-1)**(i%2)
                j += (-1)**(j%2)
                print i,j,'*'
                d = n.conjugate(d)
        if i%2==1 and j%2==0:
            d = n.conjugate(d); 
            print i,j,'**'   
        p = (uvw,t,(i,j))
        return p,d
        sys.stdout.flush()
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc,
        append2hist='PSA_WTF_CONJ: Unbunging the correlator output.')
    del(uvo)