#! /usr/bin/env python
import aipy as a, numpy as n, os, sys

aa = a.phs.ArrayLocation(('38:25:59.24',  '-79:51:02.1')) # Green Bank, WV

rewire = { 
    0 : {'x': 0, 'y':16},
    1 : {'x': 1, 'y':17},
    2 : {'x': 2, 'y':18},
    3 : {'x': 3, 'y':19},
    4 : {'x': 4, 'y':20},
    5 : {'x': 5, 'y':21},
    6 : {'x': 6, 'y':22},
    7 : {'x': 7, 'y':23},
    8 : {'x': 8, 'y':24},
    9 : {'x': 9, 'y':25},
    10: {'x':10, 'y':26},
    11: {'x':11, 'y':27},
    12: {'x':12, 'y':28},
    13: {'x':13, 'y':29},
    14: {'x':14, 'y':30},
    15: {'x':15, 'y':31},
}
bad_bls = { # 7/31/10 from Dave MacMahon
    a.miriad.ij2bl(3,12): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(3,13): ['xx'],
    a.miriad.ij2bl(4,13): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(4,14): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(4,15): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(5,14): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(5,15): ['xx','xy','yx','yy'],
    a.miriad.ij2bl(6,15): ['xx','xy','yx','yy'],
}
omit = [0,10,15]
omit += [1,2,7,8,9,16,17,18,23]
#ch0,ch1 = 600,600+1024
ch0,ch1 = 0,2048
sfreq = 0.10009765625
sdf = 97.65e-6 / 2
row1 = [ 0, 1, 2, 3, 4, 5, 6, 7]
row2 = [ 8, 9,10,11,12,13,14,15]
row3 = [16,17,18,19,20,21,22,23]
afield = [24,25,26,27,28,29,30,31]
# 0,15 noise, 10 malfunctioning
# 11_(17,18),12_(17,18,19),14_(17,18,19,20,21) != (8,9)_(17,18,19),11_19,12_(20,21,22)

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
        ## Toss bad baselines
        #bl = a.miriad.ij2bl(i,j)
        #if bad_bls.has_key(bl) and (p1+p2) in bad_bls[bl]: return p,None,None
        # Toss 4th pol combo for auto-correlations
        if i == j and (p1,p2) == ('y','x'): return p, None, None

        #t -= 2.0/24 # Recorded jultimes are 2hrs ahead?
        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
            curtime = t
        
        # Throw out chans we're not interested in right now
        d, f = d[ch0:ch1], f[ch0:ch1]
        # Correct known correlator conjugations
        if (j-i) > 8:
            d = n.conjugate(d); 
            if p1 != p2:
                p1,p2 = p2,p1
                d = n.conjugate(d)
        if (p1,p2)==('y','x'): d = n.conjugate(d)

        ni,nj = rewire[i][p1],rewire[j][p2]
        if ni in omit or nj in omit: return p, None, None
        # only save correlations between rows
        if (ni in row1 and nj in row1) or \
           (ni in row2 and nj in row2) or \
           (ni in row3 and nj in row3): return p, None, None
        if ni > nj: ni,nj,d = nj,ni,n.conjugate(d)

        return (crd,t,(ni,nj)),d,f
            
    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'sdf': sdf,
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'sfreq':sfreq + ch0*sdf,
        'freq':sfreq + ch0*sdf,
        'inttime':200000000./1024/128/4096,
        'nchan': ch1-ch0,
        'nants': 33,
        'ngains': 64,
        'nspect0': 32,
        'pol':a.miriad.str2pol['yy'],
        'telescop':'PAPER',
        #'bandpass': n.ones((33,1024), dtype=n.complex64).flatten(),
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

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='CORRECT: '+' '.join(sys.argv)+'\n'+
        'CORRECT: '+str(rewire)+'\n'+
        'CORRECT: '+str(omit)+'\n')
    del(uvo)
