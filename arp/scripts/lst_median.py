#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, ephem

o = optparse.OptionParser()
o.add_option('-b', '--bin', dest='bin', type='float', default=0.001,
    help='Bin size in LST.  Default 0.001')
o.add_option('-l', '--lstrng', dest='lstrng', default='0_6.2832',
    help='Range of LSTs to bin.')
o.add_option('--tfile', dest='tfile', type='float', default=3600,
    help='Length of time spanned by each input file.  Helps in filtering out files that are not needed for the lst range being processed.')
opts, args = o.parse_args(sys.argv[1:])

lsts = n.arange(0, 2*n.pi, opts.bin)
# Select only some of these LST bins to analyze
lstrng = map(float, opts.lstrng.split('_'))
if lstrng[0] < lstrng[1]:
    lsts = lsts.compress(n.logical_and(lsts >= lstrng[0], lsts < lstrng[1]))
else:
    lsts = lsts.compress(n.logical_or(lsts >= lstrng[0], lsts < lstrng[1]))

tfile = 2 * n.pi * opts.tfile / (24. * 3600)
dat = {}
for lst in lsts:
    dat[lst] = {}
crds = {}
jd_start = None
djd = opts.bin / (2*n.pi) * a.const.sidereal_day * ephem.second

print 'Filtering input files for LSTs of interest'
nargs = []
for f in args:
    uv = a.miriad.UV(f)
    start_t = uv['lst']
    end_t = (start_t + tfile) % (2*n.pi)
    if start_t < end_t:
        if lstrng[0] < lstrng[1]:
            if end_t < lstrng[0] or start_t > lstrng[1]: continue
        else: 
            if end_t < lstrng[0] and start_t > lstrng[1]: continue
    else:
        if lstrng[0] < lstrng[1]:
            if start_t > lstrng[1] and end_t < lstrng[0]: continue
        # Never bail if both wrap...
    nargs.append(f)

for f in nargs:
    uv = a.miriad.UV(f)
    print 'Reading', f
    sys.stdout.flush()
    if lstrng[0] < lstrng[1]:
        uv.select('lst', lstrng[0], lstrng[1])
    else:
        uv.select('lst', lstrng[1], lstrng[0], include=False)

    # Gather data from file
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        crds[bl] = uvw
        lst = n.floor(uv['lst'] / opts.bin) * opts.bin
        # Only take this LST if we have a bin for it already allocated
        if not dat.has_key(lst): continue
        try: dat[lst][bl].append(d)
        except(KeyError): dat[lst][bl] = [d]
        if jd_start is None and lst == lsts[0]:
            # Back out any fraction of an LST bin to figure a start JD
            jd_start = t - (uv['lst'] - lsts[0]) / (2*n.pi) * a.const.sidereal_day * ephem.second

uvi = a.miriad.UV(args[0])
filename = 'lst.%7.5f.uv' % jd_start
print 'Writing to', filename
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uvi)

flags = n.zeros(uvi['nchan'], dtype=n.int)
for cnt, t in enumerate(lsts):
    # Shift to center of LST bin
    ts = t + opts.bin / 2
    # Also shift JD to center of LST bin
    jd = jd_start + (cnt+.5) * djd
    lst = ephem.hours(ts)
    print 'LST:', lst, '(%f)' % ts, ' -> JD:', jd
    sys.stdout.flush()
    uvo['lst'], uvo['ra'], uvo['obsra'] = ts, ts, ts
    if len(dat[t]) == 0:
        del(dat[t])
        continue
    for bl in dat[t]:
        i,j = a.miriad.bl2ij(bl)
        preamble = (crds[bl], jd, (i,j))
        d = n.array(dat[t][bl])
        m = n.median(d, axis=0)
        m.shape = (1, m.size)
        s = n.sqrt(n.median(n.abs(d - m)**2, axis=0)).clip(1e-4, 1)
        # Should smooth s ...
        #flags = n.where(s > .1, 1, 0)
        wgts = n.exp(-2 * n.abs(d - m)**2 / s**2)
        wm = n.sum(d * wgts, axis=0) / n.sum(wgts, axis=0)
        uvo.write(preamble, wm, flags)
del(uvo)
print 'Finished writing', filename

'''
# Rebin for smoother sigma
SIGBIN = .01
nbins = int(2*n.pi/SIGBIN)
swgts,sbins = n.histogram(lsts, range=(0,2*n.pi), bins=nbins)
svals,sbins = n.histogram(lsts, weights=dat_sig, range=(0,2*n.pi), bins=nbins)
sig_bin = svals/swgts.clip(1,n.inf)

# Do weighted average of data
dat_wgt_mea = []
for i,t in enumerate(lst_t):
    ind = n.floor(t / SIGBIN).clip(0,sig_bin.size-1).astype(n.int)
'''

