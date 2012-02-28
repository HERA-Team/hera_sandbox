#! /usr/bin/env python
import aipy as a, numpy as n, os, sys, glob

def grid_jd(jds, temps, binsize=120):
    jdbin = binsize * a.ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

import optparse
o = optparse.OptionParser()
o.add_option('-t', '--tempdir', dest='tempdir',
    help='Directory containing temperature data from the (labjack) gainometer.')
opts,args = o.parse_args(sys.argv[1:])

if opts.tempdir != None:
    RECVR,ANT58,GOM_E = 0,1,2
    files = glob.glob(opts.tempdir + '/2011*.txt'); files.sort()
    lines = [L.split() for f in files for L in open(f).readlines()]
    jds = n.array([a.phs.ephem2juldate(a.ephem.date(' '.join(L[:2]))) for L in lines])
    jds -= 2*a.ephem.hour # Adjust for local time being 2 hrs ahead of UTC
    dat = n.array([map(float,L[2:]) for L in lines])
    T_r, bins = grid_jd(jds, dat[:,RECVR])
    T_c, bins = grid_jd(jds, dat[:,ANT58])
    T_b, bins = grid_jd(jds, dat[:,ANT58])
    T_l, bins = grid_jd(jds, dat[:,GOM_E])

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

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    _curtime = 0
    def _mfunc(uv, p, d, f):
        global _curtime
        crd,t,(i,j) = p
        p1,p2 = a.miriad.pol2str[uv['pol']]
        # prevent multiple entries arising from xy and yx on autocorrelations
        if i == j and (p1,p2) == ('y','x'): return p, None, None
        i,j = pol2to1[i][p1], pol2to1[j][p2]
        if i > j: i,j,d = j,i,n.conjugate(d)
        if t != _curtime:
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            _curtime = t
        
        # Deal with correlator drop-outs
        f = n.logical_or(f, n.where(n.abs(d) == 0, 1, 0))
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
            n.array([ 147.93527, 337.12580, 264.51808]) + n.array([  0.0    ,  0.0    , 0.0    ]),             
            n.array([ -122.020309402, -268.023014624, -199.815623761,]),             
            n.array([ 178.182697611, -281.228145672, 314.256932747,]),             
            n.array([ -27.26782,-368.33731, -35.49412]) + n.array([ -0.29182,  0.29581,-0.06470]),             
            n.array([ -136.271417956, -63.8538192223, -223.81696088,]),             
            n.array([ -185.394536198, 61.2146784634, -307.596987596,]),             
            n.array([  81.62552,-398.52015, 151.70308]) + n.array([ -0.90456,  0.31155,-0.03229]),             
            n.array([  61.18198, 423.06854, 115.99715]) + n.array([  0.04571,  0.00224, 0.06347]),             
            n.array([ 151.744912834, -230.476259414, 266.100117346]),             
            n.array([-121.15655,-329.37685,-197.06224]) + n.array([  0.07513,  0.00316,-0.14179]),            
            n.array([ -31.90961,-420.62509, -43.56460]) + n.array([ -0.54834,  0.39122,-0.24124]),            
            n.array([-181.51436,-188.96090,-301.06291]) + n.array([ -0.01202,  0.07766,-0.02666]),            
            n.array([ 160.72973, 208.27653, 286.82157]) + n.array([ -0.00755,  0.06142,-0.02066]),            
            n.array([ -77.8707662359, 265.443210445, -122.866728234,]),            
            n.array([  93.49461, 405.98665, 171.30307]) + n.array([  0.01422, -0.06648, 0.09750]),            
            n.array([ 137.15781,-345.85204, 247.10124]) + n.array([ -0.38124,  0.29203,-0.17590]),            
            n.array([ 72.4889166113, -363.847833856, 135.520493153]),            
            n.array([ -169.390511149, 113.82049335, -280.122814575]),            
            n.array([ -174.729494724, -51.6560694198, -289.267372556]),            
            n.array([ 33.674525626, -74.865368323, 69.2972343811,]),            
            n.array([ 221.456027623, -108.390005006, 391.362891776]),            
            n.array([ 210.58399,-180.86686, 373.06684]) + n.array([ -1.00517,  0.11158,-0.28632]),            
            n.array([ -55.3656528503, -405.662993034, -84.3890558675]),
            n.array([ -72.5558848039, 377.393814897, -113.678876716]),
            n.array([ -90.34611,   4.21680,-144.20799]) + n.array([  0.01805,  0.10158,-0.09630]),
            n.array([ -24.79049,-153.74222, -31.28959]) + n.array([ -0.35529,  0.45195,-0.14708]),
            n.array([ 210.10061544, 187.975673579, 370.683657743]),
            n.array([ -21.84807, 312.60274, -26.72816]) + n.array([  0.04485,  0.29455,-0.25724]),
            n.array([ -18.8534131847, 166.106071174, -21.0403737068]),
            n.array([  88.43837, -21.20718, 162.84732]) + n.array([  0.06337,  0.22358,-0.09150]),
            n.array([-136.73690, 313.93707,-223.87046]) + n.array([  0.00693,  0.14085,-0.14483]),
            n.array([ 230.19780727, 47.1563467525, 406.219047894,]),
            n.array([-121.89484, 147.01857,-178.53714]) + n.array([-12.34044,-21.29123,-1.62317]),
            n.array([ 121.35534,-319.42959, 209.57574]) + n.array([-13.46091,-21.08721,-2.04411]),
            n.array([  -1.18600, 298.79078,  -1.57273]) + n.array([-12.15281,-21.09783,-1.02046]),
            n.array([-150.75421,-224.46078,-258.59405]) + n.array([-13.33191,-21.31901,-2.09637]),
            n.array([-148.16634, 285.39014,-254.15270]) + n.array([-12.14134,-20.96592,-1.06129]),
            n.array([ 65.4745545234, -398.989793689, 135.937771903]),
            n.array([ 183.23862, 145.04638, 314.99738]) + n.array([-12.53417,-21.24913,-1.65605]),
            n.array([ 201.11005, 270.60894, 345.38803]) + n.array([-12.49116,-21.34547,-1.58327]),
            n.array([-187.75317, 101.63458,-322.33070]) + n.array([-13.027  ,-20.968  ,-1.872  ]),
            n.array([  32.85944,-311.36127,  57.49240]) + n.array([-13.45655,-21.32004,-2.25554]),
            n.array([ 111.79179,-360.75226, 193.12456]) + n.array([-13.38514,-21.23868,-1.85556]),
            n.array([ 185.29648,  12.47387, 318.94840]) + n.array([-12.90096,-21.24585,-1.81436]),
            n.array([  66.84088, 269.98916, 115.13990]) + n.array([-12.32095,-21.10367,-1.21838]),
            n.array([ 208.32754,-181.02402, 358.71376]) + n.array([-13.77787,-21.21314,-2.35549]),
            n.array([ 222.40198, 114.55998, 382.32980]) + n.array([-12.61714,-21.02361,-1.42653]),
            n.array([  87.62899,-157.64380, 134.49244]) + n.array([-12.57070,-20.86601,-1.67331]),
            n.array([-123.36405,   7.56840,-211.39198]) + n.array([-12.63378,-21.12748,-1.07278]),
            n.array([  42.32481,-394.59655,  73.80015]) + n.array([-13.86259,-21.23736,-2.03213]),
            n.array([ 155.42810, 103.98180, 267.54514]) + n.array([-12.63677,-21.21007,-1.64005]),
            n.array([   4.00271, 454.85825,   7.08648]) + n.array([-11.63465,-20.89803,-1.78326]),
            n.array([28.8976025965, 358.993678483, 69.7233597137,]),
            n.array([ 247.12661,  75.95094, 404.94444]) + n.array([-13.14037,-20.72110,-1.99411]),
            n.array([ 195.692326034, 148.9067559, 358.382640028]),
            n.array([  22.16270, 221.12001,  38.44946]) + n.array([-13.027  ,-20.968  ,-1.872  ]),
            n.array([ -85.96290, 360.45682,-147.01823]) + n.array([-11.24916,-21.18504,-1.59575]),
            n.array([ -22.18217, 447.51766, -37.58554]) + n.array([-11.87513,-21.21941,-1.20133]),
            n.array([ -40.13290,-349.20766, -68.17466]) + n.array([-13.60638,-21.21031,-2.56979]),
            n.array([ -38.86438, 362.86645, -66.27003]) + n.array([-11.69585,-21.02930,-1.18640]),
            n.array([ 121.811369523, 377.231448971, 229.602800651]),
            n.array([ -94.7836220969, -297.64232068, -141.49713256]),
            n.array([-161.60804, 226.51205,-277.24339]) + n.array([-12.05381,-20.93013,-1.23832]),
            n.array([170.275635,-299.76472, 293.55448]) + n.array([-13.44311,-21.21483,-2.31634]),
        ])).flatten(),
    }
    uvo.init_from_uv(uvi, override=override)
    uvo.add_var('t_recvr', 'r')
    uvo.add_var('t_cable', 'r')
    uvo.add_var('t_balun', 'r')
    uvo.add_var('t_load', 'r')
    if opts.tempdir != None:
        curtime = None
        def mfunc(uv, p, d, f):
            global curtime, uvo
            (crd, t, (i,j)) = p
            if curtime != t:
                curtime = t
                b = int(n.floor((t - bins[0]) / (bins[1] - bins[0])))
                if not n.isnan(T_r[b]): uvo['t_recvr'] = T_r[b]
                if not n.isnan(T_c[b]): uvo['t_cable'] = T_c[b]
                if not n.isnan(T_b[b]): uvo['t_balun'] = T_b[b]
                if not n.isnan(T_l[b]): uvo['t_load']  = T_l[b]
            return _mfunc(uv, p, d, f)
    else: mfunc = _mfunc

    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='CORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
