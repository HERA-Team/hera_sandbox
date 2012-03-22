#!/opt/local/Library/Frameworks/Python.framework/Versions/2.6/Resources/Python.app/Contents/MacOS/Python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse 

o = optparse.OptionParser()
o.set_usage('extract_dly_bins.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])
dly_rng = map(int, opts.dw.split('_'))

for uvfile in args:
    uvofile = uvfile + 'E'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    def mfunc(uv, p, d, f):
        crd,t,(i,j) = p
        if i == j: return p,None,None
        print t, (i,j)
        bl = a.miriad.ij2bl(i,j)
        val = n.logical_not(f).astype(n.float)
        #gain = n.sqrt(n.average(val**2))
        if not n.all(val == 0):
            # Maybe do this more than once?
            ker = n.fft.ifft(val)
            _d = n.fft.ifft(d * val)
            _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            _d = info['res']
            #_d[opts.dw+1:-opts.dw] = 0
            #d -= n.fft.fft(_d) * val
            _d[:dly_rng[0]] = 0
            _d[dly_rng[1]:-dly_rng[1]+1] = 0
            _d[-dly_rng[0]+1:] = 0
            d = n.fft.fft(_d) * val
        return p, d, f

    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    # Apply the pipe to the data
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='FILTER_SKY: dw=%d clean=%f\n' % (opts.dw, opts.clean))

