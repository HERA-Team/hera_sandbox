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
    return [date] + temps

def grid_jd(jds, temps, binsize=120):
    jdbin = binsize * ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

not_conj_bls = [
    a.miriad.ij2bl(0,6),
    a.miriad.ij2bl(0,7),
    a.miriad.ij2bl(1,6),
    a.miriad.ij2bl(1,7),
    a.miriad.ij2bl(2,4),
    a.miriad.ij2bl(2,5),
    a.miriad.ij2bl(3,4),
    a.miriad.ij2bl(3,5),
]

def mfunc(uv, p, d, f):
    crd,t,(i,j) = p
    bl = a.miriad.ij2bl(i,j)
    if bl in not_conj_bls: return p, d, f
    else: return p, n.conj(d), f

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    override = {
        'sfreq':0.100390625,
        'freq':0.100390625,
        'inttime':0.747520029545,
        }
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='Conjugate most baselines.\nRewrote sfreq/freq/inttime to be correct\n')
    del(uvo)

