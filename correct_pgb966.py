#! /usr/bin/env python

import aipy as a, numpy as n, sys, os, ephem, optparse
from time import strftime, gmtime, mktime, strptime


def tstr2jd(tstr, ifmt='%m/%d/%y %H:%M:%S', tz='UTC'):
    try:
        tstr = strftime('%Y/%m/%d %H:%M:%S',
            gmtime(mktime(strptime(tstr+' '+tz, ifmt+' %Z'))))
        return a.phs.ephem2juldate(ephem.date(tstr))
    except(ValueError): return []

def parse_gom_line(line):
    fields = line.split()
    date = tstr2jd(' '.join(fields[:2]))
    temps = map(float, fields[2:])
    return [date] + temps

def grid_jd(jds, temps, binsize=600):
    jdbin = binsize * ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

offset = 630
nchan = 1024
chans = numpy.arange(offset, offset+nchan)
skip_ants = [7]
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
    aipy.miriad.ij2bl(0,5),
    aipy.miriad.ij2bl(0,6),
    aipy.miriad.ij2bl(1,6),
    aipy.miriad.ij2bl(0,7),
    aipy.miriad.ij2bl(1,7),
    aipy.miriad.ij2bl(2,7),
]

def mfunc(uv, p, d, f):
    uvw, t, (i,j) = p
    bl = aipy.miriad.ij2bl(i,j)
    p1,p2 = aipy.miriad.pol2str[uv['pol']]
    if p1 != p2 and bl in swap_xy_yx_pol_labels:
        p1,p2 = p2,p1
        d = numpy.conjugate(d)
    ni = rewire[i][p1]
    nj = rewire[j][p2]
    if i == j and (p1,p2) == ('y','x'): return p, None, None
    if ni in skip_ants or nj in skip_ants: return p, None, None
    # Fix X-Engine conj bug
    if j - i > 4: d = numpy.conjugate(d)
    # Fix 75% data loss in X-Engine loopback bug
    if j - i > 4: d /= 0.75
    if ni > nj:
        ni,nj = nj,ni
        if (p1,p2) != ('y','x'): d = numpy.conjugate(d)
    elif ni < nj and (p1,p2) == ('y','x'):
        d = numpy.conjugate(d)
    p = (uvw,t,(ni,nj))
    return p, d.take(chans), f.take(chans)

o = optparse.OptionParser()
o.add_option('-t', '--tempdir', dest='tempdir',
    help='Directory containing temperature data from the gainometer.')
opts,args = o.parse_args(sys.argv[1:])

if opts.tempdir != None:
    dat = []
    for f in glob.glob(opts.tempdir + '/Temp*.txt'):
        lines = [parse_gom_line(L) for L in open(f).readlines()]
        dat += [L for L in lines if len(L) == 5]
    dat = n.array(dat)
    T_r, bins = grid_jd(dat[:,0], dat[:,1])
    T_c, bins = grid_jd(dat[:,0], dat[:,2])
    T_b, bins = grid_jd(dat[:,0], dat[:,3])
    T_l, bins = grid_jd(dat[:,0], dat[:,4])
        

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(filename+'c', status='new')
    uvo.add_var('t_recvr', 'r')
    uvo.add_var('t_cable', 'r')
    uvo.add_var('t_balun', 'r')
    uvo.add_var('t_load', 'r')
    sfreq = uvi['sfreq'] + offset * uvi['sdf']
    bp = uvi['bandpass']
    bp.shape = (2048, 8)
    bp = numpy.concatenate([bp, bp], axis=1)
    bp = bp.take(chans, axis=0)
    bp = bp.transpose()
    override = {'pol':aipy.miriad.str2pol['yy'],
        'nchan':nchan, 'nchan0':nchan, 'sfreq':sfreq,
        'freq':sfreq, 'bandpass':bp.flatten(), 'nants':16, 'npol':1,
        'freqs':(16,) + (nchan, sfreq, uvi['sdf']) * 16}
    uvo.init_from_uv(uvi, override=override)
    if opts.tempdir != None:
        curtime = None
        def mfunc(uv, p, d, f):
            global curtime
            (crd, t, (i,j)) = p
            if curtime != t:
                curtime = t
                b = int(n.floor((t - bins[0]) / (bins[1] - bins[0])))
                uv['t_recvr'] = T_r[b]
                uv['t_cable'] = T_c[b]
                uv['t_balun'] = T_b[b]
                uv['t_load']  = T_l[b]
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='Relabeled as 16 ants, conjugated, trimmed freqs, rm ant#7\n')
    del(uvo)

