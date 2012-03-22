#! /usr/bin/env python
import aipy as a, numpy as n, os, sys

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m

pol2to1 = {
    0: {'x': 0, 'y': 1},
    1: {'x': 2, 'y': 3},
    2: {'x': 4, 'y': 5},
    3: {'x': 6, 'y': 7},
    4: {'x': 8, 'y': 9},
    5: {'x': 10, 'y': 11},
    6: {'x': 12, 'y': 13},
    7: {'x': 14, 'y': 15},
    8: {'x': 16, 'y': 17},
    9: {'x': 18, 'y': 19},
    10: {'x': 20, 'y': 21},
    11: {'x': 22, 'y': 23},
    12: {'x': 24, 'y': 25},
    13: {'x': 26, 'y': 27},
    14: {'x': 28, 'y': 29},
    15: {'x': 30, 'y': 31},
    16: {'x': 32, 'y': 33},
    17: {'x': 34, 'y': 35},
    18: {'x': 36, 'y': 37},
    19: {'x': 38, 'y': 39},
    20: {'x': 40, 'y': 41},
    21: {'x': 42, 'y': 43},
    22: {'x': 44, 'y': 45},
    23: {'x': 46, 'y': 47},
    24: {'x': 48, 'y': 49},
    25: {'x': 50, 'y': 51},
    26: {'x': 52, 'y': 53},
    27: {'x': 54, 'y': 55},
    28: {'x': 56, 'y': 57},
    29: {'x': 58, 'y': 59},
    30: {'x': 60, 'y': 61},
    31: {'x': 62, 'y': 63},
}

for filename in sys.argv[1:]:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0
    def mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        p1,p2 = a.miriad.pol2str[uv['pol']]
        # prevent multiple entries arising from xy and yx on autocorrelations
        if i == j and (p1,p2) == ('y','x'): return p, None, None

        i,j = pol2to1[i][p1], pol2to1[j][p2]

        if i > j: i,j,d = j,i,n.conjugate(d)

        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            curtime = t
        
        return (crd,t,(i,j)),d,f

    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        #'sdf': sdf,
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        #'sfreq': sfreq,
        #'freq': sfreq,
        #'inttime': 5.37,
        'nchan': 1024,
        'nants': 64,
        'ngains': 128,
        'nspect0': 64,
        'pol':a.miriad.str2pol['xx'],
        'telescop':'PAPER',
        #'bandpass': n.ones((33,1024), dtype=n.complex64).flatten(),
        'antpos': n.transpose(n.array([
            [0., 0., 0.], #0
            [0., 0., 0.], #1 
            [0., 0., 0.], #2
            [0., 0., 0.], #3
            [0., 0., 0.], #4
            [0., 0., 0.], #5
            [0., 0., 0.], #6
            [0., 0., 0.], #7
            [0., 0., 0.], #8
            [0., 0., 0.], #9
            [0., 0., 0.], #10
            [0., 0., 0.], #11
            [0., 0., 0.], #12
            [0., 0., 0.], #13
            [0., 0., 0.], #14
            [0., 0., 0.], #15
            [0., 0., 0.], #16
            [0., 0., 0.], #17
            [0., 0., 0.], #18
            [0., 0., 0.], #19
            [0., 0., 0.], #20
            [0., 0., 0.], #21
            [0., 0., 0.], #22
            [0., 0., 0.], #23
            [0., 0., 0.], #24
            [0., 0., 0.], #25
            [0., 0., 0.], #26
            [0., 0., 0.], #27
            [0., 0., 0.], #28
            [0., 0., 0.], #29
            [0., 0., 0.], #30
            [0., 0., 0.], #31
            [0., 0., 0.], #32
            [0., 0., 0.], #33
            [0., 0., 0.], #34
            [0., 0., 0.], #35
            [0., 0., 0.], #36
            [0., 0., 0.], #37
            [0., 0., 0.], #38
            [0., 0., 0.], #39
            [0., 0., 0.], #40
            [0., 0., 0.], #41
            [0., 0., 0.], #42
            [0., 0., 0.], #43
            [0., 0., 0.], #44
            [0., 0., 0.], #45
            [0., 0., 0.], #46
            [0., 0., 0.], #47
            [0., 0., 0.], #48
            [0., 0., 0.], #49
            [0., 0., 0.], #50
            [0., 0., 0.], #51
            [0., 0., 0.], #52
            [0., 0., 0.], #53
            [0., 0., 0.], #54
            [0., 0., 0.], #55
            [0., 0., 0.], #56
            [0., 0., 0.], #57
            [0., 0., 0.], #58
            [0., 0., 0.], #59
            [0., 0., 0.], #60
            [0., 0., 0.], #61
            [0., 0., 0.], #62
            [0., 0., 0.], #63
        ])).flatten(),
    }

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='CORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
