#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, ephem
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--lst_res', type='float', default=10.,
    help='Resolution in seconds for binning in LST.  Default is 10.')
o.add_option('--lst_rng', default='0_6.2832',
    help='Range of LSTs to bin.')
o.add_option('--flag_thresh', type='float', default=0.5,
    help='Fraction of data that must be unflagged a bin for the output to be unflagged.')
o.add_option('--tfile', type='float', default=600,
    help='Length of time spanned by each input file.  Helps in filtering out files that are not needed for the lst range being processed.')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

# Select only some of these LST bins to analyze
opts.lst_rng = map(float, opts.lst_rng.split('_'))
def in_lst_range(lst):
    if opts.lst_rng[0] < opts.lst_rng[1]:
        return lst >= opts.lst_rng[0] and lst < opts.lst_rng[1]
    else:
        return lst >= opts.lst_rng[0] or lst < opts.lst_rng[1]

def lstbin(lst):
    lst_res = opts.lst_res / a.const.sidereal_day * (2*n.pi)
    return C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst,lst_res=lst_res),lst_res=lst_res)[-1]

lstbins = n.arange(0, 2*n.pi, 2*n.pi*opts.lst_res/a.const.sidereal_day)
lstbins = [lstbin(lst) for lst in lstbins if in_lst_range(lst)]

dat,cnt = {}, {}
for lst in lstbins: dat[lst],cnt[lst] = {}, {}
crds = {}
jd_start = None

print 'Filtering input files for LSTs of interest'
nargs = []
for f in args:
    uv = a.miriad.UV(f)
    (crd,t,bl),_d,_f = uv.read(raw=True)
    aa.set_jultime(t)
    start_t = aa.sidereal_time()
    aa.set_jultime(t + opts.tfile * a.ephem.second)
    end_t = aa.sidereal_time()
    if start_t < end_t:
        if opts.lst_rng[0] < opts.lst_rng[1]:
            if end_t < opts.lst_rng[0] or start_t > opts.lst_rng[1]: continue
        else:
            if end_t < opts.lst_rng[0] and start_t > opts.lst_rng[1]: continue
    else:
        if opts.lst_rng[0] < opts.lst_rng[1]:
            if start_t > opts.lst_rng[1] and end_t < opts.lst_rng[0]: continue
        # Never bail if both wrap...
    nargs.append(f)

jds = {}
for f in nargs:
    uv = a.miriad.UV(f)
    print 'Reading', f
    sys.stdout.flush()
    # Unsure whether to trust lsts of uv files
    #if lstrng[0] < lstrng[1]:
    #    uv.select('lst', lstrng[0], lstrng[1])
    #else:
    #    uv.select('lst', lstrng[1], lstrng[0], include=False)

    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Gather data from file
    curtime = None
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            aa.set_jultime(t)
            lst = lstbin(aa.sidereal_time())
            if dat.has_key(lst): jds[lst] = min(jds.get(lst,n.Inf), t)
        # Only take this LST if we have a bin for it already allocated
        if not dat.has_key(lst): continue
        bl = a.miriad.ij2bl(i,j)
        crds[bl] = uvw
        dat[lst][bl] = dat[lst].get(bl,0) + n.where(f,0,d)
        cnt[lst][bl] = cnt[lst].get(bl,0) + n.logical_not(f).astype(n.int)

lsts = [lst for lst in dat if len(dat[lst]) > 0]
lsts.sort()
lst_start = lsts[0]
jd_start = jds[lst_start]
djd_dlst = a.const.sidereal_day / (2*n.pi) * a.ephem.second

uvi = a.miriad.UV(args[0])
filename = 'lst.%7.5f.uv' % jd_start
print 'Writing to', filename
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uvi)

for lst in lsts:
    t = jd_start + (lst - lst_start) * djd_dlst
    print 'LST:', a.ephem.hours(lst), '(%f)' % lst, ' -> JD:', t
    sys.stdout.flush()
    uvo['lst'], uvo['ra'], uvo['obsra'] = lst, lst, lst
    for bl in dat[lst]:
        i,j = a.miriad.bl2ij(bl)
        preamble = (crds[bl], t, (i,j))
        cmax = n.max(cnt[lst][bl])
        d = dat[lst][bl] / cnt[lst][bl].clip(1, n.Inf)
        f = n.where(cnt[lst][bl] < cmax * opts.flag_thresh, 1, 0)
        uvo.write(preamble, d, f)
del(uvo)
print 'Finished writing', filename

