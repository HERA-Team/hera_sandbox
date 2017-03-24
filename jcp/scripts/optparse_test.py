#! /usr/bin/env python

import aipy as a, numpy as n
import capo as C
import optparse, sys, os
import glob

POL_WGTS = {
    'I': {'xx': 1. , 'yy': 1. },
    'Q': {'xx': 1. , 'yy':-1. },
    'U': {'xy': 1. , 'yx': 1. },
    'V': {'xy':-1.j, 'yx': 1.j},
}

LIN2STOKES = {
    'xx':'I',
    'yy':'Q',
    'xy':'U',
    'yx':'V',
}

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--stokes', help='Comma-delimited list of Stokes parameters to form (e.g. I,Q,U,V).')
o.add_option('--xx', default=None, help='Filenames for the xx pols.')
o.add_option('--yy', default=None, help='Filenames for the yy pols.')
o.add_option('--xy', default=None, help='Filenames for the xy pols.')
o.add_option('--yx', default=None, help='Filenames for the yx pols.')
opts,args = o.parse_args(sys.argv[1:])

stokes = opts.stokes.split(',')
pol = {}
for s in stokes: pol.update(POL_WGTS[s])
pol = ','.join(pol.keys())

if opts.xx == None and opts.yy == None and opts.xy == None and opts.yx == None:
    files = args
    NPOLFILES = 1
elif opts.xx != None and opts.yy != None and opts.xy == None and opts.yx == None:
    if opts.stokes != 'I':
        print 'Invalid Stokes parameter selection for the input linear pols'
        sys.exit()
    else:
        xxfiles = glob.glob(opts.xx)
        yyfiles = glob.glob(opts.yy)
        files = zip(xxfiles,yyfiles)
        NPOLFILES = 2
elif opts.xx != None and opts.yy != None and opts.xy != None and opts.yx != None:
    xxfiles = glob.glob(opts.xx)
    yyfiles = glob.glob(opts.yy)
    xyfiles = glob.glob(opts.xy)
    yxfiles = glob.glob(opts.yx)
    files = zip(xxfiles,yyfiles,xyfiles,yyfiles)
    NPOLFILES = 4

for filelist in files:
    if NPOLFILES > 1: outfile = os.path.commonprefix(filelist) + 'P'
    else: outfile = filelist + 'P'
    print filelist, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dsum,dwgt = {}, {}
    for s in stokes: dsum[s], dwgt[s] = {}, {}
    curtime = None
    print '    Converting %s to %s' % (pol, ','.join(stokes))
    if NPOLFILES > 1:
        uvlist = []
        for linpol in xrange(NPOLFILES):
            uvi = a.miriad.UV(filelist[linpol])
            ants = a.scripting.parse_ants(opts.ant, uvi['nants'])
            a.scripting.uv_selector(uvi, ants)
            for (crd,t,(i,j)),d,f in uvi.all(raw=True):
                bl = a.miriad.ij2bl(i,j)
                if t != curtime:
                    curtime = t
                    #aa.set_jultime(t)
                    for s in stokes:
                        dsum[s][t], dwgt[s][t] = {}, {}
                for s in stokes:
                    try: wgt = POL_WGTS[s][a.miriad.pol2str[uvi['pol']]]
                    except(KeyError): continue
                    dsum[s][t][bl] = dsum[s][t].get(bl, 0) + n.where(f, 0, wgt*d)
                    dwgt[s][t][bl] = dwgt[s][t].get(bl, 0) + n.abs(wgt)*n.logical_not(f).astype(n.int)
    
print dsum
