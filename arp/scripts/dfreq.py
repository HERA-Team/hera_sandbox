#!/opt/local/Library/Frameworks/Python.framework/Versions/2.6/Resources/Python.app/Contents/MacOS/Python
"""
Rotate zenith UV data to a particular source.  Can specify 'zen' to phase data
to zenith, or nothing at all to just remove delay/offset phase components.
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('dfreq.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--width', dest='width', type='int', default=25,
    help='Width in channels of boxcar filter to use for removing smooth freq components.')
opts,args = o.parse_args(sys.argv[1:])

_filter = None
def filter(d, f, width=opts.width):
    global _filter
    if _filter is None:
        dim = d.size
        _filter = n.zeros_like(d)
        #_filter[dim/2-width/2:dim/2+width/2+1] = 1
        _filter[:width/2+1] = 1
        _filter[-width/2:] = 1
        _filter /= _filter.sum()
        _filter = n.fft.fft(_filter)
    d = n.where(f, 0, d)
    _f = f.astype(n.float)
    _f = n.fft.ifft(_filter * n.fft.fft(_f))
    d -= n.fft.ifft(_filter * n.fft.fft(d)) / _f.clip(1,n.Inf)
    #d = n.fft.ifft(_filter * n.fft.fft(d)) / f.clip(1,n.Inf)
    f = n.logical_or(f, n.where(_f >= .25, 1, 0))
    d = n.where(f > 0, 0, d)
    return d, f

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
if not opts.src is None:
    if not opts.src.startswith('zen'):
        srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
        src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]
    else: src = 'z'
else: src = None
del(uv)

# A pipe to use for phasing to a source
curtime = None
def phs(uv, p, d, f):
    global curtime
    uvw, t, (i,j) = p
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        if not src is None and not type(src) == str: src.compute(aa)
    if i == j: return p, d, f
    try:
        phs = aa.gen_phs(src, i, j, resolve_src=False)
        d *= phs
        d,f = filter(d,f)
        d /= phs
    except(a.phs.PointingError): d *= 0
    return p, d, f

# Process data
for filename in args:
    uvofile = filename + 'F'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs, raw=True)

