#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('-n', '--nint', dest='nint', type='int',
    help='Number of integrations to add together into output file.')
opts,args = o.parse_args(sys.argv[1:])

#assert(opts.nint > 0)
i_rec = 0
dat, cnt, tbuf = {}, {}, {}
def mfunc(uv, p, d, f):
    global dat, cnt, tbuf
    global n_rec,i_rec
    crd,t,(i,j) = p
    bl = a.miriad.ij2bl(i,j)
    if not dat.has_key(bl):
        dat[bl] = 0; cnt[bl] = 0; tbuf[bl] = []
    dat[bl] += n.where(f,0,d)
    cnt[bl] += n.logical_not(f).astype(n.int)
    tbuf[bl].append(t)
    print i_rec
    if i_rec == n_rec-1:
        f = n.where(cnt[bl] > n_rec/2, 0, 1)
        d = n.where(f, 0, dat[bl]/cnt[bl])
        t = n.average(tbuf[bl])
        p = crd,t,(i,j)
        print len(d),n.max(d)
        del(dat[bl])
    else:
        d, f = None, None
        i_rec += 1
    return p, d, f

for filename in args:
    print filename, '->', filename+'A'
    if os.path.exists(filename+'A'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    n_times = 0
    t_min = 0
    t_max = 0
    n_rec = 0
    for (uvw,t,(i,j)),d in uvi.all():
        if t_min ==0 or t<t_min: t_min = t
        if t_max==0 and t>t_max: t_max = t
        n_rec +=1
    uvi.rewind()
    uvo = a.miriad.UV(filename+'A', status='new')
    override = {
        'inttime':t_max-t_min,
    }
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='SUM_INT: Added together all %d integrations over %d seconds. '%(n_rec,t_max-t_min))
    del(uvo)
    
