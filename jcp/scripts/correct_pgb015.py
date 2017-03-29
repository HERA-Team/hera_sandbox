#! /usr/bin/env python

import aipy as a, numpy as n, sys, os, ephem, optparse, glob
from time import strftime, gmtime, mktime, strptime


def tstr2jd(tstr, ifmt='%m/%d/%y %H:%M:%S', tz=''):
    try:
        #tstr = strftime('%Y/%m/%d %H:%M:%S',
        #    gmtime(mktime(strptime(tstr+' '+tz, ifmt+' %Z'))))
        tstr = strftime('%Y/%m/%d %H:%M:%S', strptime(tstr, ifmt))
        return a.phs.ephem2juldate(ephem.date(tstr))
    except(ValueError): return []

def parse_gom_line(line, filename):
    fields = line.split()
    if len(fields) < 6 and len(fields) > 0:
        ifmt='%m%d%y %H:%M:%S'
        date = tstr2jd(' '.join([filename[-10:-4], fields[0]]), ifmt=ifmt)
        temps = map(float, fields[1:])
    else:
        date = tstr2jd(' '.join(fields[:2]))
        temps = map(float, fields[2:])
    print date
    return [date] + temps

def grid_jd(jds, temps, binsize=120):
    jdbin = binsize * ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

offset = 630
nchan = 1024
chans = n.arange(offset, offset+nchan)
#skip_ants = [7, 13, 15]
skip_ants = []
rewire = {
    6: {'x': 8, 'y': 9},
    7: {'x':10, 'y':11},
    4: {'x': 4, 'y': 5},
    5: {'x': 6, 'y': 7},
    2: {'x': 0, 'y': 1},
    3: {'x': 2, 'y': 3},
    0: {'x':12, 'y':13},
    1: {'x':14, 'y':15},
}

swap_xy_yx_pol_labels = [
    a.miriad.ij2bl(0,5),
    a.miriad.ij2bl(0,6),
    a.miriad.ij2bl(1,6),
    a.miriad.ij2bl(0,7),
    a.miriad.ij2bl(1,7),
    a.miriad.ij2bl(2,7),
]

def _mfunc(uv, p, d, f):
    uvw, t, (i,j) = p
    bl = a.miriad.ij2bl(i,j)
    p1,p2 = a.miriad.pol2str[uv['pol']]
    if p1 != p2 and bl in swap_xy_yx_pol_labels:
        p1,p2 = p2,p1
        d = n.conjugate(d)
    ni = rewire[i][p1]
    nj = rewire[j][p2]
    if i == j and (p1,p2) == ('y','x'): return p, None, None
    if ni in skip_ants or nj in skip_ants: return p, None, None
    # Fix X-Engine conj bug
    if j - i > 4: d = n.conjugate(d)
    # Fix 75% data loss in X-Engine loopback bug
    if j - i > 4: d /= 0.75
    if ni > nj:
        ni,nj = nj,ni
        if (p1,p2) != ('y','x'): d = n.conjugate(d)
    elif ni < nj and (p1,p2) == ('y','x'):
        d = n.conjugate(d)
    p = (uvw,t,(ni,nj))
    return p, d.take(chans), f.take(chans)

o = optparse.OptionParser()
o.add_option('-t', '--tempdir', dest='tempdir',
    help='Directory containing temperature data from the gainometer.')
opts,args = o.parse_args(sys.argv[1:])

if opts.tempdir != None:
    dat = []
    files = glob.glob(opts.tempdir + '/Temp*.txt')
    files.sort()
    for f in files:
        print 'Reading', f
        lines = [parse_gom_line(L, f) for L in open(f).readlines()]
        dat += [L for L in lines if len(L) == 5]
    dat = n.array(dat)
    T_r, bins = grid_jd(dat[:,0], dat[:,1])
    T_c, bins = grid_jd(dat[:,0], dat[:,2])
    T_b, bins = grid_jd(dat[:,0], dat[:,3])
    T_l, bins = grid_jd(dat[:,0], dat[:,4])
    #import pylab as p
    #p.plot(bins[:-1], T_r)
    #p.show()
        

histline = ' '.join([s for s in sys.argv if not s in args])
for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    sfreq = uvi['sfreq'] + offset * uvi['sdf']
    bp = uvi['bandpass']
    bp.shape = (2048, 8)
    bp = n.concatenate([bp, bp], axis=1)
    bp = bp.take(chans, axis=0)
    bp = bp.transpose()
    override = {'pol':a.miriad.str2pol['yy'],
        'nchan':nchan, 'nchan0':nchan, 'sfreq':sfreq,
        'freq':sfreq, 'bandpass':bp.flatten(), 'nants':16, 'npol':1,
        'freqs':(16,) + (nchan, sfreq, uvi['sdf']) * 16}
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
                b = int(n.floor((t - bins[0]) / (bins[1] - bins[0])))
                if curtime is None and n.isnan(T_r[b]):
                    # If this is the first integration in file, find the closest bin that works
                    b_off = 1
                    while n.isnan(T_r[b-b_off]) and n.isnan(T_r[b+b_off]):
                        b_off += 1
                    if not n.isnan(T_r[b-b_off]): b -= b_off
                    else: b += b_off
                if not n.isnan(T_r[b]): uvo['t_recvr'] = T_r[b]
                if not n.isnan(T_c[b]): uvo['t_cable'] = T_c[b]
                if not n.isnan(T_b[b]): uvo['t_balun'] = T_b[b]
                if not n.isnan(T_l[b]): uvo['t_load']  = T_l[b]
                curtime = t
            return _mfunc(uv, p, d, f)
    else: mfunc = _mfunc
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=histline+' %s\n' % filename)
    del(uvo)

