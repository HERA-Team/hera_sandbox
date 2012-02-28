#! /usr/bin/env python

import aipy as a, numpy as n, os, optparse, sys

aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9'))



rewire_run1 = { 
    0: {'x':0, 'y':1},
    1: {'x':2, 'y':3},
    2: {'x':4, 'y':5},
    3: {'x':6, 'y':7},
    4: {'x':8, 'y':9},
    5: {'x':10, 'y':11},
    6: {'x':12, 'y':13},
    7: {'x':14, 'y':15},
    8: {'x':16, 'y':17},
    9: {'x':18, 'y':19},
    10: {'x':20, 'y':21},
    11: {'x':22, 'y':23},
    12: {'x':24, 'y':25},
    13: {'x':26, 'y':27},
    14: {'x':28, 'y':29},
    15: {'x':30, 'y':31},
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
        p1,p2 = a.miriad.pol2str[uv['pol']]
        #if i == j and (p1,p2) == ('y','x'): return p, None, None

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
            
            if (p1,p2) != ('y','x'): d = n.conjugate(d)
        elif ni < nj and (p1,p2) == ('y','x'):
            d = n.conjugate(d)

        p = crd,t,(ni,nj)
#        print t,i,j,a.miriad.pol2str[uv['pol']],'->',ni,nj,'xx'
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
        #'bandpass': n.resize(uvi['bandpass'], (64,1024)).flatten(),
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
#    if is_run1(uvi['time']): override['pol'] = a.miriad.str2pol['yy']

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc,
        append2hist='Rewrote sfreq/freq/inttime to be correct?\nSet initial antenna positions based on PSA-32 data.')
    del(uvo)
