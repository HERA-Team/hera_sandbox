#! /usr/bin/env python
import aipy as a, numpy as n,os
import sys, optparse, ephem
import capo as C

'''
This is an update to the lst binner code that now produces meta data in the uv
file, such as the number of integrations that went into each bin.  This allows
the code to properly add new data to previously lst binned files.
'''

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, src=True)
o.add_option('--lst_res', type='float', default=10.,
    help='Resolution in seconds for binning in LST.  Default is 10.')
o.add_option('--lst_rng', default='0_23.999',
    help='Range of LSTs to bin, hours. [default=0_23.999]')
o.add_option('--tfile', type='float', default=600,
    help='Length of time spanned by each input file.  Helps in filtering out files that are not needed for the lst range being processed.')
o.add_option('--altmax', type='float', default=0,
    help="Maximum allowed altitude of source, in degrees, before data are omitted.  Handy for omitting Sun data.  Default is 0.")
o.add_option('--stats', default='all',
    help="Statistics to include in meta data.  Options are 'all' (default), 'none', 'cnt' (the number of integrations in each lst bin), 'min' and/or 'max' (the min and max amplitude of the visibilities in each lst bin', 'median' (the median value in each lst bin), and 'var' (the variance in each lst bin).  Multiple values can be chosen with commas e.g. '--stats=min,max,var'.")
o.add_option('--median', action='store_true', dest='median', default=False,
    help="Use a median filter to remove outliers from each lst bin.")
o.add_option('--nsig', type='float', default=3.,
    help="Number of sigma outlier to flag in median filter.")
opts, args = o.parse_args(sys.argv[1:])

# ---- Functions for lst binning ----

def in_lst_range(lst):
# Returns True/False if a given lst in within the range specifed by opts.lstrng
    if opts.lst_rng[0] < opts.lst_rng[1]:
        return lst >= opts.lst_rng[0] and lst < opts.lst_rng[1]
    else:
        return lst >= opts.lst_rng[0] or lst < opts.lst_rng[1]

def lstbin(lst):
# Chooses an lst bin for a given lst
    lst_res = opts.lst_res / a.const.sidereal_day * (2*n.pi)
    return C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst,lst_res=lst_res),lst_res=lst_res)[-1]

# ---- Code ----

# Reads the first file to create an antenna array
uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

# Creates a src object if opts.src flag is used
src = None
if not opts.src is None:
    srclist,cutoff,catalog = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalog)
    src = cat.values()[0]

# Parses the stats flag
opts.stats = map(str, opts.stats.split(','))

# Create the lst bins inside the desired range
opts.lst_rng = map(lambda x: n.pi/12*(float(x)), opts.lst_rng.split('_')) #Hours
print opts.lst_rng
lstbins = n.arange(0, 2*n.pi, 2*n.pi*opts.lst_res/a.const.sidereal_day)
lstbins = [lstbin(lst) for lst in lstbins if in_lst_range(lst)]

# Create dat dict with subdicts for each lst bin
dat = {}
for lst in lstbins: dat[lst] = {}
crds = {}
jd_start = None

# Checks each input file for lsts in range, returns nargs = list of files to use
print 'Filtering input files for LSTs of interest'
nargs = []
for f in args:
    uv = a.miriad.UV(f)
    # Read the first timestamp from the file
    (crd,t,bl),_d,_f = uv.read(raw=True)
    aa.set_jultime(t)
    if not src is None:
        src.compute(aa)
        src_alt_start = src.alt
    start_t = aa.sidereal_time()
    # Set the time to the end of the file (file length specified by opts.tfile)
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
        # ARP: Never bail if both wrap...
        # JCP: I have yet to wrap my head around how/why this would happen
    # Include the file is src is below altmax at beginning or end
    if src is None or (src_alt_start < opts.altmax or src_alt_end < opts.altmax):
        nargs.append(f)

# Places the data into lst bins, but does not actually combine or average.
jds = {}
files = {}
for filename in nargs:
    uv = a.miriad.UV(filename)
    print 'Reading', filename
    # Determine if file is lst binned already
    if 'cnt' in uv.vars(): 
        lstbinned = True
        print 'Input file is lst binned'
    else: 
        lstbinned = False
        print 'Input file is from one observation'
    sys.stdout.flush()
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Gather data from file
    curtime = None
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            aa.set_jultime(t)
            if not src is None: src.compute(aa)
            lst = lstbin(aa.sidereal_time())
            if dat.has_key(lst):
                jds[lst] = jds.get(lst,[]) + [t] # Keep track of all jds that contribute
                files[lst] = files.get(lst,[]) + [filename] # Keep track of all files that contribute
            curtime = t
        # Only include this integration if it falls within the defined range
        if not dat.has_key(lst): continue
        # Don't include the integration is src is above altmax
        if not src is None and src.alt >= opts.altmax: continue
        # Records coords of baseline in dict, for writing into uv file later
        blp = a.pol.ijp2blp(i,j,uv['pol'])
        crds[blp] = uvw
        # Keep track of how many (unflagged) samples go into the data
        # If input file already lstbinned, weight integrations and place in bin
        # Otherwise, place in bin without weights
        if not dat[lst].has_key(blp): dat[lst][blp] = {}
        if lstbinned:
            dat[lst][blp]['cnt'] = dat[lst][blp].get('cnt',[]) + [uv['cnt']]
            dat[lst][blp]['vis'] = dat[lst][blp].get('vis',[]) + [uv['cnt']*n.where(f,0,d)]
        else:
            dat[lst][blp]['cnt'] = dat[lst][blp].get('cnt',[]) + [n.logical_not(f)]
            dat[lst][blp]['vis'] = dat[lst][blp].get('vis',[]) + [n.where(f,0,d)]

# Check that data actually got written
lsts = [lst for lst in dat if len(dat[lst]) > 0] # only record bins with data

lsts.sort()
if len(lsts) == 0:
    print 'No LST bins with data.  Exiting...'
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
jd_start = jd_start + (lst_start - lsts[0]) * djd_dlst
lst_start = lsts[0]

# Initialize the output file
#XXX this section of code assumes that all the input files look like the first one.  this will deliver BAD functionality if input files are a mix of lst binned and non lst binned files
uvi = a.miriad.UV(args[0])
invars = uvi.vars() #remember what variables were in the input files
filename=os.path.basename(args[0])
# DCJ: This is for beamformed .bm_<srcname> files
if filename.split('.')[-1].startswith('bm'):
    filename='lst.%7.5f.uv.%s' % (jd_start,filename.split('.')[-1])
else: filename = 'lst.%7.5f.uv' % jd_start
print 'Writing to', filename
if os.path.exists(filename):
    print filename,"exists"
    sys.exit(1)
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uvi)

# What to write in the case of zero data in an lst bin
dzero = n.zeros(uvi['nchan'], dtype=n.complex64)
fzero = n.ones(uvi['nchan'], dtype=n.int)

# Add the variables for the statistics if needed
if opts.stats == ['all']: opts.stats = ['cnt','min','max','median','var']
if opts.stats == ['none']: opts.stats = []
if 'cnt' in opts.stats and 'cnt' not in uvo.vars(): uvo.add_var('cnt', 'd')
if 'min' in opts.stats and 'min' not in uvo.vars(): uvo.add_var('min', 'd')
if 'max' in opts.stats and 'max' not in uvo.vars(): uvo.add_var('max', 'd')
if 'median' in opts.stats and 'median' not in uvo.vars(): uvo.add_var('median', 'd')
if 'var' in opts.stats and 'var' not in uvo.vars(): uvo.add_var('var', 'd')

# This section does the binning and statistics 
for lst in lsts:
    t = jd_start + (lst - lst_start) * djd_dlst
    print 'LST:', a.ephem.hours(lst), '(%f)' % lst, ' -> JD:', t
    sys.stdout.flush()
    uvo['lst'], uvo['ra'], uvo['obsra'] = lst, lst, lst
    for blp in blps:
        i,j,uvo['pol'] = a.pol.blp2ijp(blp)
        preamble = (crds[blp], t, (i,j))
        try:
            d = n.array(dat[lst][blp]['vis'])
            w = n.array(dat[lst][blp]['cnt'])
            d = n.ma.array(d, mask=n.where(d==0, 1, 0))
            # Applies median filter
            if opts.median:
                d_med = n.ma.median(d, axis=0)
                d_med.shape = (1,) + d_med.shape
                d_res = n.ma.abs(d - d_med)
                d_sig = n.ma.median(d_res, axis=0)
                d_sig.shape = (1,) + d_sig.shape
                d = n.ma.masked_where(d_res > opts.nsig * d_sig, d)
            # Calculate the statistics as needed
            if 'cnt' in opts.stats: cnt = n.sum(w, axis=0)
            if 'min' in opts.stats: dmin = n.ma.min(n.abs(d), axis=0).filled(0)
            if 'max' in opts.stats: dmax = n.ma.max(n.abs(d), axis=0).filled(0)
            if 'median' in opts.stats: median = n.ma.median(n.abs(d), axis=0).filled(0)
            if 'var' in opts.stats: var = n.ma.var(d, axis=0).filled(0)
            # Do the averaging
            #d = n.ma.mean(d, axis=0).filled(0)
            #Weighted average
            d = n.ma.sum(w*d,axis=0)/n.sum(w,axis=0)
            f = n.where(d == 0, 1, 0)
       # This happens if we are missing data for a desired LST bin
        except(KeyError):
            d,f = dzero, fzero
            if 'cnt' in opts.stats: cnt = n.zeros_like(fzero)
            if 'minmax' in opts.stats: dmin, dmax = n.zeros_like(fzero), n.zeros_like(fzero)
            if 'median' in opts.stats: median = n.zeros_like(fzero)
            if 'var' in opts.stats: var = n.zeros_like(fzero)
        # Set the statistics variables in the uv object and write the data
        # Commands if input files DO NOT have vars to begin with
        if 'cnt' in opts.stats and 'cnt' not in invars: uvo['cnt'] = cnt.astype(n.double)
        if 'min' in opts.stats and 'min' not in invars: uvo['min'] = dmin.astype(n.double)
        if 'max' in opts.stats and 'max' not in invars: uvo['max'] = dmax.astype(n.double)
        if 'median' in opts.stats and 'median' not in invars: uvo['median'] = median.astype(n.double)
        if 'var' in opts.stats and 'var' not in invars: uvo['var'] = var.astype(n.double)
        # Commands if input files ALREADY have vars
        if 'cnt' in opts.stats and 'cnt' in invars: uvo['cnt'] += cnt.astype(n.double)
        if 'min' in opts.stats and 'min' in invars: uvo['min'] = n.minimum(uvo['min'],dmin.astype(n.double))
        if 'max' in opts.stats and 'max' in invars: uvo['max'] = n.maximum(uvo['max'],dmax.astype(n.double))
        #if 'median' in opts.stats and 'median' in invars: uvo['median'] += median.astype(n.double) #XXX no idea what to do here
        #if 'var' in opts.stats and 'var' in invars: uvo['var'] += var.astype(n.double) #XXX this one should be doable?

        uvo.write(preamble, d, f)

del(uvo)
print 'Finished writing', filename

