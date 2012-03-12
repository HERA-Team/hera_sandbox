#!/usr/bin/env python
#
#  correct_poco335_psapol.py
#  
#
#  Created by Danny Jacobs on 5/19/10.
#  PAPER Project
#

import aipy as a, numpy as n, os, optparse, sys

aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9'))

def is_run1(jd): return (jd > 2455229.00 and jd < 2455233.00)

#rewire_run1 = { 
#    0: {'x':0, 'y':1},
#    1: {'x':2, 'y':3},
#    2: {'x':4, 'y':5},
#    3: {'x':6, 'y':7},
#    4: {'x':8, 'y':9},
#    5: {'x':10, 'y':11},
#    6: {'x':12, 'y':13},
#    7: {'x':14, 'y':15},
#    8: {'x':16, 'y':17},
#    9: {'x':18, 'y':19},
#    10: {'x':20, 'y':21},
#    11: {'x':22, 'y':23},
#    12: {'x':24, 'y':25},
#    13: {'x':26, 'y':27},
#    14: {'x':28, 'y':29},
#    15: {'x':30, 'y':31},
#}
rewire_ant = {
6:3,
7:3,
0:2,
1:2,
2:1,
3:1,
4:0,
5:0
}
rewire_pol = {
6:'x',
7:'y',
0:'x',
1:'y',
2:'x',
3:'y',
4:'x',
5:'y'
}

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0
    def mfunc(uv, p, d):
        crd,t,(i,j) = p
        ni = rewire_ant[i]
        nj = rewire_ant[j]
        pi,pj = rewire_pol[i],rewire_pol[j]

        if t!= curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()        

#        else: return p, None
        if ni > nj:
            ni,nj = nj,ni
            pi,pj = pj,pi
            if (pi,pj) != ('y','x'): d = n.conjugate(d)
        elif ni < nj and (pi,pj) == ('y','x'):
            d = n.conjugate(d)
        uvo['pol'] = a.miriad.str2pol[pi+pj]
        p = crd,t,(ni,nj)
#        print i,j,'xx','->',ni,nj,pi,pj
        return p,d
            
    override = {
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'sfreq':0.10009765625,
        'freq':0.10009765625,
        'inttime':200000000./1024/128/1024,
        'nants': 4,
        'ngains': 8,
        'nspect0': 4,
        #'bandpass': n.resize(uvi['bandpass'], (64,1024)).flatten(),
        'antpos': n.transpose(n.array([
            [58.053568,423.509161,115.997158],
        ])).flatten(),
    }
#    if is_run1(uvi['time']): override['pol'] = a.miriad.str2pol['yy']

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc,
        append2hist='Rewrote sfreq/freq/inttime to be correct?\nSet initial antenna positions based on PSA-32 data.')
    del(uvo)
