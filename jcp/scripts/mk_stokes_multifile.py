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

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--stokes', help='Comma-delimited list of Stokes parameters to form (e.g. I,Q,U,V).')
o.add_option('--xx', default=None, help='Filenames for the xx pols.')
o.add_option('--yy', default=None, help='Filenames for the yy pols.')
o.add_option('--xy', default=None, help='Filenames for the xy pols.')
o.add_option('--yx', default=None, help='Filenames for the yx pols.')
o.add_option('--mode', default='paper', help='Makes assumptions about the filename format.  Options are "paper" and "fhd".')
opts,args = o.parse_args(sys.argv[1:])

if not opts.mode == 'paper' and not opts.mode == 'fhd':
    print 'Invalid mode selected... Exiting.'
    sys.exit()

stokes = opts.stokes.split(',')
pol = {}
for s in stokes: pol.update(POL_WGTS[s])
pol = ','.join(pol.keys())
#uv = a.miriad.UV(args[0])
#ants = a.scripting.parse_ants(opts.ant, uv['nants'])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#del(uv)

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
    files = zip(xxfiles,yyfiles,xyfiles,yxfiles)
    NPOLFILES = 4
else:
    print 'Invalid combination of linear pols... Exiting.'
    sys.exit()

print NPOLFILES

for filelist in files:
    if NPOLFILES > 1: 
        if opts.mode == 'fhd': outfile = os.path.commonprefix(filelist) + 'P'
        if opts.mode == 'paper': outfile = os.path.commonprefix(filelist) + os.path.basename(filelist[0]).split('.')[-1] + 'P'
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
        polcnt = 0
        for linpol in xrange(NPOLFILES):
            polcnt += 1
            uvi = a.miriad.UV(filelist[linpol])
            ants = a.scripting.parse_ants(opts.ant, uvi['nants'])
            a.scripting.uv_selector(uvi, ants)
            for (crd,t,(i,j)),d,f in uvi.all(raw=True):
                bl = a.miriad.ij2bl(i,j)
                if t != curtime:
                    curtime = t
                    #aa.set_jultime(t)
                    if polcnt == 1:
                        for s in stokes:
                            dsum[s][t], dwgt[s][t] = {}, {}
                for s in stokes:
                    try: wgt = POL_WGTS[s][a.miriad.pol2str[uvi['pol']]]
                    except(KeyError): continue
                    dsum[s][t][bl] = dsum[s][t].get(bl, 0) + n.where(f, 0, wgt*d)
                    dwgt[s][t][bl] = dwgt[s][t].get(bl, 0) + n.abs(wgt)*n.logical_not(f).astype(n.int)
            uvi.rewind()
    else:
        uvi = a.miriad.UV(filelist)
        #a.scripting.uv_selector(uvi, ants, pol)
        ants = a.scripting.parse_ants(opts.ant, uvi['nants'])
        a.scripting.uv_selector(uvi, ants)
        #uvi.select('decimate', opts.decimate, opts.decphs)
        for (crd,t,(i,j)),d,f in uvi.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if t != curtime:
                curtime = t
                #aa.set_jultime(t)
                for s in stokes:
                    dsum[s][t], dwgt[s][t] = {}, {}
            for s in stokes:
                print s, POL_WGTS[s][a.miriad.pol2str[uvi['pol']]]
                try: wgt = POL_WGTS[s][a.miriad.pol2str[uvi['pol']]]
                except(KeyError): continue
                dsum[s][t][bl] = dsum[s][t].get(bl, 0) + n.where(f, 0, wgt*d)
                dwgt[s][t][bl] = dwgt[s][t].get(bl, 0) + n.abs(wgt)*n.logical_not(f).astype(n.int)
        uvi.rewind()

    print '    Writing output file'
    uvo = a.pol.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'pol':a.miriad.str2pol['I']})
    for pol in dsum.keys():
        uvo.write_pol(pol)
        times = n.array(dsum[pol].keys()).astype(float)
        times.sort()
        for t in times:
            for bl in dsum[pol][t].keys():
                i,j = a.miriad.bl2ij(bl)
                p = (n.array([0.,0.,0.]),t,(i,j))
                #try: 
                _dsum,_dwgt = dsum[pol][t].pop(bl), dwgt[pol][t].pop(bl)
                f = n.where(_dwgt <= 1., 1, 0) #flag things that are flagged in any linear pol
                wgt = _dwgt.clip(1,n.Inf)
                d = _dsum / wgt
                #except(KeyError): 
                #    d, f = None, None
                uvo.write(p,d,f)
            
    uvo._wrhd('history',uvo['history'] + 'COMBINE_POL:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




