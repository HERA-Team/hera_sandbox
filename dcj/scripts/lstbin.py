#! /usr/bin/env python
import aipy as a, numpy as n,os
import sys, optparse, ephem
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, src=True)
o.add_option('--lst_res', type='float', default=10.,
    help='Resolution in seconds for binning in LST.  Default is 10.')
o.add_option('--lst_rng', default='0_6.2832',
    help='Range of LSTs to bin.')
o.add_option('--flag_thresh', type='float', default=0.5,
    help='Fraction of data that must be unflagged a bin for the output to be unflagged.')
o.add_option('--tfile', type='float', default=600,
    help='Length of time spanned by each input file.  Helps in filtering out files that are not needed for the lst range being processed.')
o.add_option('--altmax', type='float', default=0,
    help="Maximum allowed altitude of source, in degrees, before data are omitted.  Handy for omitting Sun data.  Default is 0.")
o.add_option('--nogaps', action='store_true',
    help='Record a spectrum for every LST in the chosen range, even if a blank one must be created to fill a gap.')

opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)
src = None
if not opts.src is None:
    srclist,cutoff,catalog = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalog)
    src = cat.values()[0]

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
def lst2date(lst,date):
    "Output the next julian date corresponding to input lst forward from input date."
    obj = a.fit.RadioFixedBody(lst,0)
    obj.compute(aa)
    aa.set_jultime(n.floor(date)) #start the day at noon UTC
    aa.update()
    aa.date =  aa.next_transit(obj)
    aa.update()
    return aa.get_jultime()
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
    if not src is None:
        src.compute(aa)
        src_alt_start = src.alt
    start_t = aa.sidereal_time()
    aa.set_jultime(t + opts.tfile * a.ephem.second)
    if not src is None:
        src.compute(aa)
        src_alt_end = src.alt
    start_t = aa.sidereal_time()
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
    if src is None or (src_alt_start < opts.altmax or src_alt_end < opts.altmax):
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
            if not src is None: src.compute(aa); print src.alt
            lst = lstbin(aa.sidereal_time())
            if dat.has_key(lst): jds[lst] = min(jds.get(lst,n.Inf), t)
            curtime = t
        # Only take this LST if we have a bin for it already allocated
        if not dat.has_key(lst): continue
        if not src is None and src.alt >= opts.altmax: continue
        blp = a.pol.ijp2blp(i,j,uv['pol'])
        crds[blp] = uvw
        dat[lst][blp] = dat[lst].get(blp,0) + n.where(f,0,d)
        cnt[lst][blp] = cnt[lst].get(blp,0) + n.logical_not(f).astype(n.int)

if opts.nogaps:
    lsts = lstbins # record all bins
else:
    lsts = [lst for lst in dat if len(dat[lst]) > 0] # only record bins with data
lsts.sort()
if len(lsts) == 0:
    print 'No LST bins with data.  Exitting...'
    sys.exit(0)
# Get list of all blps to record, just to make sure that each time has a record for each blp
blps = {}
for lst in dat:
    for blp in dat[lst]:
        blps[blp] = None
blps = blps.keys()
# Find a starting jd for recording in the file
lst_start, jd_start = n.Inf, n.Inf
for lst in jds:
    if jds[lst] < jd_start:
        lst_start, jd_start = lst, jds[lst]
djd_dlst = a.const.sidereal_day / (2*n.pi) * a.ephem.second
jd_start = jd_start + (lst_start - lsts[0]) * djd_dlst
lst_start = lsts[0]
#lst_start = lsts[0]
#jd_start = jds[lst_start]
uvi = a.miriad.UV(args[0])
filename=os.path.basename(args[0])
if filename.find('bm')>0:
    filename='lst.%7.5f.uv.%s' % (jd_start,filename.split('.')[-1])
else:filename = 'lst.%7.5f.uv' % jd_start
print 'Writing to', filename
if os.path.exists(filename):
    print filename,"exists"
    sys.exit(1)
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uvi)
# XXX could think about adding a variable that keeps track of how many integrations went into a bin

dzero = n.zeros(uvi['nchan'], dtype=n.complex64)
fzero = n.ones(uvi['nchan'], dtype=n.int)
for lst in lsts:
    #t = jd_start + (lst - lst_start) * djd_dlst
    t = lst2date(lst,jd_start)
    print 'LST:', a.ephem.hours(lst), '(%f)' % lst, ' -> JD:', t
    sys.stdout.flush()
    uvo['lst'], uvo['ra'], uvo['obsra'] = lst, lst, lst
    #for blp in dat[lst]:
    for blp in blps:
        i,j,uvo['pol'] = a.pol.blp2ijp(blp)
        preamble = (crds[blp], t, (i,j))
        try:
            cmax = n.max(cnt[lst][blp])
            d = dat[lst][blp] / cnt[lst][blp].clip(1, n.Inf)
            f = n.where(cnt[lst][blp] < cmax * opts.flag_thresh, 1, 0)
        except(KeyError): # This happens if we are missing data for a desired LST bin
            d,f = dzero, fzero
        uvo.write(preamble, d, f)
del(uvo)
print 'Finished writing', filename

