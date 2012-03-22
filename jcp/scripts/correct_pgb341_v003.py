#! /usr/bin/env python

'''A correct scripts for pgb341.  Once again, this has been thrown together by JP so is likely terrible.'''

import aipy as a, numpy as n, os, optparse, sys

aa = a.phs.ArrayLocation(('38:25:59.24',  '-79:51:02.1')) # Green Bank, WV


pol2to1 = { 
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

no_ants = [7]

ant2corrinp = {
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    5: 5,
    6: 6,
    8: 7,
    9: 8,
    10: 9,
    11: 10,
    12: 11,
    13: 12,
    14: 13,
    15: 14,
    16: 15,
    17: 16,
    18: 17,
    19: 18,
    20: 19,
    21: 20,
    22: 21,
    23: 22,
    24: 23,
    25: 24,
    26: 25,
    27: 26,
    28: 27,
    29: 28,
    30: 29,
    31: 30,
    32: 31,
}

ant2corr = {
    0: 0,
    1: 2,
    2: 4,
    3: 6,
    4: 8,
    5: 10,
    6: 12,
    8: 14,
    9: 16,
    10: 18,
    11: 20,
    12: 22,
    13: 24,
    14: 26,
    15: 28,
    16: 30,
    17: 1,
    18: 3,
    19: 5,
    20: 7,
    21: 9,
    22: 11,
    23: 13,
    24: 15,
    25: 17,
    26: 19,
    27: 21,
    28: 23,
    29: 25,
    30: 27,
    31: 29,
    32: 31,
}

corr2ant = {
    0: 0,
    1: 17,
    2: 1,
    3: 18,
    4: 2,
    5: 19,
    6: 3,
    7: 20,
    8: 4,
    9: 21,
    10: 5,
    11: 22,
    12: 6,
    13: 23,
    14: 8,
    15: 24,
    16: 9,
    17: 25,
    18: 10,
    19: 26,
    20: 11,
    21: 27,
    22: 12,
    23: 28,
    24: 13,
    25: 29,
    26: 14,
    27: 30,
    28: 15,
    29: 31,
    30: 16,
    31: 32
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
        global curtime
        crd,t,(i,j) = p
        t += -2.0/24 #timezone error?
        p1,p2 = a.miriad.pol2str[uv['pol']]
        if i == j and (p1,p2) == ('y','x'): return p, None
        #i = pol2to1[i][p1]
        #j = pol2to1[j][p2]
        if t!= curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
            curtime = t
        
        #i = ant2corr[i]
        #j = ant2corr[j]
        if (j-i)>8:
            d = n.conjugate(d); 
            if p1 != p2:
                p1,p2 = p2,p1
                d = n.conjugate(d)
        if (p1,p2)==('y','x'): d = n.conjugate(d)
        
        i = pol2to1[i][p1]
        j = pol2to1[j][p2]
        if i > j:
            i,j = j,i
            d = n.conjugate(d)


        #i = corr2ant[i]
        #j = corr2ant[j]
        #if i > j:
        #    i, j = j, i
        #    d = n.conjugate(d)
        p = crd,t,(i,j)
        #d[orbchan].mask = 1
        #uv['bandpass'][no_ants,:] == 0j
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
        'nants': 33,
        'ngains': 64,
        'nspect0': 32,
        'pol':a.miriad.str2pol['xx'],
        'telescop':'PAPER',
        'bandpass': n.resize(uvi['bandpass'], (33,1024)).flatten(),
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
        append2hist='Rewrote sdf/sfreq/freq/inttime to be correct\nSet initial antenna positions to ZERO.\nChanged times by -2 hours.\nFixed? correlator conjugation errors.')
    del(uvo)
