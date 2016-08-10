#! /usr/bin/env python
import aipy as a, numpy as n,os
import sys, optparse, ephem
import capo as C

MEDIAN = True
CUT = 0.95
NSIG = 3

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
files = {}
for filename in nargs:
    uv = a.miriad.UV(filename)
    print 'Reading', filename
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
            if not src is None: src.compute(aa)#; print src.alt
            lst = lstbin(aa.sidereal_time())
            #if dat.has_key(lst): jds[lst] = min(jds.get(lst,n.Inf), t)
            if dat.has_key(lst):
                jds[lst] = jds.get(lst,[]) + [t] # keep track of all jds that contribute
                files[lst] = files.get(lst,[]) + [filename] # keep track of all files that contribute
            curtime = t
        # Only take this LST if we have a bin for it already allocated
        if not dat.has_key(lst): continue
        if not src is None and src.alt >= opts.altmax: continue
        blp = a.pol.ijp2blp(i,j,uv['pol'])
        crds[blp] = uvw
        if MEDIAN:
            dat[lst][blp] = dat[lst].get(blp,[]) + [n.where(f,0,d)]
        else:
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
for lst in jds.keys():
    jds[lst] = n.array(jds[lst])
    files[lst] = n.array(files[lst])
    jd_min = n.min(jds[lst])
    if jd_min < jd_start:
        lst_start, jd_start = lst, jd_min
djd_dlst = a.const.sidereal_day / (2*n.pi) * a.ephem.second
jd_start = jd_start + (lsts[0] - lst_start) * djd_dlst
lst_start = lsts[0]
#lst_start = lsts[0]
#jd_start = jds[lst_start]

uvi = a.miriad.UV(args[0])
filename=os.path.basename(args[0])
if filename.split('.')[-1].startswith('bm'):
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
niter = 0.
for lst in lsts:
    bad_files = {}
    t = jd_start + (lst - lst_start) * djd_dlst
    print 'LST:', a.ephem.hours(lst), '(%f)' % lst, ' -> JD:', t
    sys.stdout.flush()
    uvo['lst'], uvo['ra'], uvo['obsra'] = lst, lst, lst
    #for blp in dat[lst]:
    for blp in blps:
        i,j,uvo['pol'] = a.pol.blp2ijp(blp)
        preamble = (crds[blp], t, (i,j))
        try:
            cmax = n.max(cnt[lst][blp])
            if MEDIAN:
                d = n.array(dat[lst][blp])
                c = n.array(cnt[lst][blp])
                if True:
                    #mask = n.where(d == 0, 1, 0)
                    #print n.sum(mask, axis=0)
                    d = n.ma.array(d, mask=n.where(d==0, 1, 0))
                    #import pylab
                    #for d_ in d: pylab.plot(d_.real, 'r.')
                    d_med = n.ma.median(d, axis=0); d_med.shape = (1,) + d_med.shape
                    #print d_med.dtype
                    #pylab.plot(d_med[0].real, 'k-')
                    d_res = n.ma.abs(d - d_med)
                    d_sig = n.ma.median(d_res, axis=0); d_sig.shape = (1,) + d_sig.shape
                    #print d_res.dtype, d_sig.dtype, d.dtype
                    #pylab.plot(d_med[0].real+NSIG*d_sig[0],'k:')
                    #pylab.plot(d_med[0].real-NSIG*d_sig[0],'k:')
                    #mask = n.where(d_res > NSIG * d_sig, 1, 0)
                    #print n.ma.sum(mask, axis=0)
                    #print mask.shape
                    d = n.ma.masked_where(d_res > NSIG * d_sig, d)
                    #print d.dtype
                    #d = n.ma.average(d, axis=0).filled(0)
                    #for d_ in d: pylab.plot(d_.real, 'b.')
                    d = n.ma.mean(d, axis=0).filled(0)
                    #print d.dtype
                    #pylab.plot(d.real, 'g')
                    #pylab.show()
                    # XXX technically might want to update cnt based on this new flagging.  not critical
                    #d = n.median(dat[lst][blp], axis=0)
                else:
                    ad = n.argsort(n.abs(d), axis=0)
                    cut = max(1, int(CUT*d.shape[0]))
                    print d.shape[0]
                    for i in xrange(d.shape[1]):
                        d[:,i] = d[:,i][ad[:,i]]
                        jd_sort = jds[lst][ad[:,i]]
                        files_sort = files[lst][ad[:,i]]
                        for f in files_sort[cut:]:
                            f = str(f)
                            bad_files[f] = bad_files.get(f,0) + 1
                    #import pylab, capo
                    #capo.arp.waterfall(d, mode='real')
                    #pylab.show()
                    d = n.sum(d[:cut],axis=0) / c.clip(1,n.Inf) / CUT # XXX hack for now to renormalize correctly
            else:
                d = dat[lst][blp] / cnt[lst][blp].clip(1, n.Inf)
            f = n.where(cnt[lst][blp] < cmax * opts.flag_thresh, 1, 0)
        except(KeyError): # This happens if we are missing data for a desired LST bin
            d,f = dzero, fzero
        niter += 1.
        uvo.write(preamble, d, f)
    bfiles = n.array([f for f in bad_files])
    bcnt = n.array([bad_files[f] for f in bfiles])
    #n.savez('bad_%f.npz' % lst, files=bfiles, cnt=bcnt)
del(uvo)
print 'Finished writing', filename
print niter
