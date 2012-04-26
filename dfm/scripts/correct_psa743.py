#! /usr/bin/env python
import aipy as a, numpy as n, os, sys

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m

cabling = {
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

mistakes = {
    40 : {'x': (40,'y'), 'y': (40,'x')}, 
    55 : {'x': (55,'y'), 'y': (55,'x')}
}

badants = [11]
badbls = [a.miriad.ij2bl(58,63),a.miriad.ij2bl(58,61),a.miriad.ij2bl(56,59)]

apply_bp = 12250000.

rfi_chans  = '0_36,49_60,78_83,86_90,102_110,114'
rfi_chans += ',121_134,168_171,187_188,341,648_649,759_776,1144'
rfi_chans += ',1540,1700_1708,1795,1826_1829,1865_1871'

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
        pi,pj = a.miriad.pol2str[uv['pol']]
        ## Toss bad baselines
        bl = a.miriad.ij2bl(i,j)
        
        #t -= 2.0/24 # Recorded jultimes are 2hrs ahead?
       	#Perform the rewire
        i,j = cabling[i],cabling[j]
        if i in mistakes.keys():
            ii = i
            i = mistakes[ii][pi][0]
            pi = mistakes[ii][pi][1]
        if j in mistakes.keys():
            jj = j
            j = mistakes[jj][pj][0]
            pj = mistakes[jj][pj][1]
        uvo['pol'] = a.miriad.str2pol[pi+pj]

        if i > j: 
            i,j = j,i
            d = n.conjugate(d)
        
        if i == j: return None,None,None
        if i in badants: return None,None,None
        if j in badants: return None,None,None
        if a.miriad.ij2bl(i,j) in badbls: return None,None,None

        d *= apply_bp
        chans = a.scripting.parse_chans(rfi_chans,uv['nchan'])
        f[chans] = 1

        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
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
        'nchan': 2048,
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
