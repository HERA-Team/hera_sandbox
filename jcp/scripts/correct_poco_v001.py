#! /usr/bin/env python

import aipy as a, numpy as n, os, optparse, sys

aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9'))

conj_bls=[]

def is_run1(jd): return (jd > 2455112.00 and jd < 2455120.00)
def is_run2(jd): return (jd > 2455228.00 and jd < 2455233.00)

rewire_run1={0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7}
rewire_run2={0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7}

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = None
    def mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        if bl in conj_bls: d = n.conj(d)
        if is_run1(t):
            #if i == 2 or j == 2: return p, None, None
            d = n.conj(d)
            i,j = rewire_run1[i], rewire_run1[j]
            uvo['pol'] = a.miriad.str2pol['xx']
        elif is_run2(t):
            #if i == 2 or j == 2: return p, None, None
	    d = n.conj(d)
	    i,j = rewire_run2[i], rewire_run2[j]
            uvo['pol'] = a.miriad.str2pol['xx']
        else: return p, None, None
        if i > j: i,j,d = j,i,n.conj(d)
        if curtime != t:
            #if is_run1(t): print 'Processing as run 1'
            #elif is_run2(t): print 'Processing as run 2'
            curtime = t
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
        p = crd,t,(i,j)
        return p,d,f
            
    override = {
	'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'sfreq':0.100390625,
        'freq':0.100390625,
        'inttime':0.747520029545,
        'nants': 32,
        'ngains': 32,
        'nspect0': 32,
        'bandpass': n.resize(uvi['bandpass'], (32,256)).flatten(),
        'antpos': n.transpose(n.array([
            [145.439661,338.209928,264.518083],
            [-121.502757,-270.051873,-201.079975],
            [175.483874,-282.046474,309.593714],
            [-24.544834,-368.528784,-35.494120],
            [-135.977107,-65.373043,-223.715356],
            [-185.770595,60.445251,-307.412232],
            [84.568610,-397.906007,151.703088],
            [58.053568,423.509161,115.997158],
            [148.405177,-231.475974,263.305593],
            [-118.718939,-330.263286,-197.062245],
            [-28.800063,-420.849441,-43.564604],
            [-180.112865,-190.297251,-301.062917],
            [159.186049,209.458742,286.821573],
            [-79.830818,266.416356,-122.828467],
            [90.491568,406.666552,171.303074],
            [139.710135,-344.828915,247.101248],
            [75.008275,-366.663944,135.807286],
            [-170.082246,113.392564,-280.090332],
            [-174.288927,-52.870923,-289.036411],
            [34.236471,-75.376892,69.315436],
            [222.470464,-108.873368,391.637799],
            [211.914965,-179.305581,373.066846],
            [-52.446982,-408.925812,-83.683533],
            [-75.327786,379.129646,-113.829018],
            [-90.374808,3.548977,-144.207995],
            [-23.653561,-153.921245,-31.289596],
            [208.418197,189.287085,370.725255],
            [-24.157809,312.432739,-26.728163],
            [-19.920924,166.712893,-21.010944],
            [88.592695,-20.552987,162.847328],
            [-139.053365,312.917932,-223.870462],
            [229.945829,48.161862,406.414507],
        ])).flatten(),
    }
    if is_run1(uvi['time']): override['pol'] = a.miriad.str2pol['xx']
    if is_run2(uvi['time']): override['pol'] = a.miriad.str2pol['xx']
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='Conjugate most baselines.\nRewrote sfreq/freq/inttime to be correct\n')
    del(uvo)
    if curtime is None:
        print 'Empty output file.  Deleting...'
        os.system('rm -rf %s' % (filename+'c'))

