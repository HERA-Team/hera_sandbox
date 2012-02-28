#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('-n', '--nint', dest='nint', type='int',
    help='Number of integrations to add together into output file.')
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

for filename in args:
    print filename, '->', filename+'I'
    if os.path.exists(filename+'I'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'I', status='new')
    override = {
        'inttime':uvi['inttime'] * opts.nint,
    }
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='SUM_INT: Added together integrations. nint=%d\n' % (opts.nint))
    del(uvo)
    
