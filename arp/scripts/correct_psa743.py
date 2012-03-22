#! /usr/bin/env python
import aipy as a, numpy as n, os, sys

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m

rewire = {
	0  : 60 ,
	1  : 61 ,
	2  : 62 ,
	3  : 63 ,
	4  : 56 ,
	5  : 57 , 
	6  : 58 ,
	7  : 59 ,
	8  : 52 ,
	9  : 53 ,
	10 : 54 ,
	11 : 55 , 
	12 : 48 ,
	13 : 49 ,
	14 : 50 ,
	15 : 51 ,
    16 : 44 ,
    17 : 45 ,
    18 : 46 ,
    19 : 47 ,
    20 : 40 ,
    21 : 41 ,
    22 : 42 ,
    23 : 43 ,
    24 : 36 ,
    25 : 37 ,
    26 : 38 ,
    27 : 39 ,
    28 : 11 ,
    29 : 33 ,
    30 : 34 ,
    31 : 35
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
        ## Toss bad baselines
        bl = a.miriad.ij2bl(i,j)
        
        #t -= 2.0/24 # Recorded jultimes are 2hrs ahead?
       	i,j = rewire[i],rewire[j]
        
        if i > j: 
            i,j = j,i
            d = n.conjugate(d)

	if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
            curtime = t
        
        # Throw out chans we're not interested in right now
        #d, f = d[ch0:ch1], f[ch0:ch1]
        
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
        'nspect0': 32,
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
        ])).flatten(),
    }

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='CORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
