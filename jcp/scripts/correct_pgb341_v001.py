#! /usr/bin/env python

'''A correct scripts for pgb341.  Once again, this has been thrown together by JP so is likely terrible.'''

import aipy as a, numpy as n, os, optparse, sys

aa = a.phs.ArrayLocation(('38:25:59.24',  '-79:51:02.1')) # Green Bank, WV


rewire_run1 = { 
    0: {'x':0, 'y':17},
    1: {'x':1, 'y':18},
    2: {'x':2, 'y':19},
    3: {'x':3, 'y':20},
    4: {'x':4, 'y':21},
    5: {'x':5, 'y':22},
    6: {'x':6, 'y':23},
    7: {'x':8, 'y':24},
    8: {'x':9, 'y':25},
    9: {'x':10, 'y':26},
    10: {'x':11, 'y':27},
    11: {'x':12, 'y':28},
    12: {'x':13, 'y':29},
    13: {'x':14, 'y':30},
    14: {'x':15, 'y':31},
    15: {'x':16, 'y':32},
}

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
orbchan = n.arange(278,287)


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
        t += -2.0/24 #timezone error?
        p1,p2 = a.miriad.pol2str[uv['pol']]

#        if is_run1(t):
        ni = rewire_run1[i][p1]
        nj = rewire_run1[j][p2]
        if t!= curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
#        else: return p, None
        if ni > nj:
            ni,nj = nj,ni
            
        #if (p1,p2) != ('y','x'): d = n.conjugate(d)
        #elif ni < nj and (p1,p2) == ('y','x'):
        d = n.conjugate(d)

        p = crd,t,(ni,nj)
        d[orbchan].mask = 1
        return p,d
            
    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'sdf': 97.65e-6,
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'sfreq':0.10009765625,
        'freq':0.10009765625,
        'inttime':200000000./1024/128/4096,
        'nants': 32,
        'ngains': 64,
        'nspect0': 32,
        'pol':a.miriad.str2pol['xx'],
        'telescop':'PAPER',
        'bandpass': n.resize(uvi['bandpass'], (64,1024)).flatten(),
        'antpos': n.transpose(n.array([
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
        ])).flatten(),
    }
#    if is_run1(uvi['time']): override['pol'] = a.miriad.str2pol['yy']

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc,
        append2hist='Rewrote sdf/sfreq/freq/inttime to be correct\nSet initial antenna positions to ZERO.\nChanged times by -2 hours.\n')
    del(uvo)
