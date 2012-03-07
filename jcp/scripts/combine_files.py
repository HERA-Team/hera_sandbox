#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os, math

o = optparse.OptionParser()
o.add_option('-n', '--nint', dest='nint', type='int',
    help='Number of integrations to add together file.')
o.add_option('--nfiles', dest='nfiles', type='int',
    help='Number of files to combine together.')
opts,args = o.parse_args(sys.argv[1:])

assert(opts.nint > 0)

dat, cnt, tbuf = {}, {}, {}
def mfunc(uv, p, d, f):
    global dat, cnt, tbuf
    crd,t,(i,j) = p
    bl = a.miriad.ij2bl(i,j)
    if not dat.has_key(bl):
        dat[bl] = 0; cnt[bl] = 0; tbuf[bl] = []
    dat[bl] += n.where(f,0,d)
    cnt[bl] += n.logical_not(f).astype(n.int)
    tbuf[bl].append(t)
    if len(tbuf[bl]) == opts.nint:
        f = n.where(cnt[bl] > opts.nint/2, 0, 1)
        d = n.where(f, 0, dat[bl]/cnt[bl])
        t = n.average(tbuf[bl])
        p = crd,t,(i,j)
        del(dat[bl])
    else:
        d, f = None, None
    return p, d, f

print len(args)

filelists = n.array_split(range(len(args)),(len(args)/opts.nfiles))

for filelist in filelists:
    print args[filelist[0]], '->', args[filelist[0]]+'M'
    if os.path.exists(args[filelist[0]]+'M'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(args[filelist[0]]+'M', status='new')
    uvtemp = a.miriad.UV(args[filelist[0]])
    override = {
        'inttime':uvtemp['inttime'] * opts.nint,
    }
    uvo.init_from_uv(uvtemp, override=override)
    del(uvtemp)

    for file in filelist:
        uvi = a.miriad.UV(args[file])
        uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='SUM_INT: Added together integrations. nint=%d\nSummed %d files together\n' % (opts.nint, opts.nfiles))

    del(uvo)
    
