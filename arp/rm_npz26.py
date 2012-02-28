#! /usr/bin/env python
import numpy as n, aipy as a
import optparse, sys, os, glob

EXP_NOISE = 1.

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-D','--DIR', dest='dir', default='.',
    help="Directory where npz files are found.")
opts,args = o.parse_args(sys.argv[1:])

try:
    import fftw3
    print 'Using FFTW FFT'
    _fftw3_dat, _fftw3_fwd, _fftw3_rev = None, None, None
    def fft2(d):
        global _fftw3_fwd, _fftw3_dat
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_fwd()
        return _fftw3_dat
    def ifft2(d):
        global _fftw3_rev, _fftw3_dat
        if _fftw3_rev is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_rev = fftw3.Plan(_fftw3_dat, None, direction='backward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_rev()
        return _fftw3_dat
except(ImportError):
    print 'Using numpy FFT'
    fft2, ifft2 = n.fft.fft2, n.fft.ifft2


def gen_filter(shape, dw, drw, ratio=.25):
    filter = n.ones(shape)
    x1,x2 = drw, -drw
    if x2 == 0: x2 = shape[0]
    y1,y2 = dw, -dw
    if y2 == 0: y2 = shape[1]
    filter[x1+1:x2,0] = 0
    filter[0,y1+1:y2] = 0
    filter[1:,1:] = 0
    x,y = n.indices(shape).astype(n.float)
    x -= shape[0]/2
    y -= shape[1]/2
    r2 = (x/(ratio*drw+.5))**2 + (y/(ratio*dw+.5))**2
    r2 = a.img.recenter(r2, (shape[0]/2, shape[1]/2))
    filter += n.where(r2 <= 1, 1, 0)
    return filter.clip(0,1)

for filename in args:
    uvi = a.miriad.UV(filename)
    uvofile = filename + 'd'
    print filename, '->', uvofile
    if os.path.exists(uvofile):
        print '    File exists: skipping'
        continue
    aa = a.cal.get_aa(opts.cal, uvi['sdf'], uvi['sfreq'], uvi['nchan'])
    srcs = {}
    srcdata = {}
    for npzfilename in glob.glob('%s/%s*.npz' % (opts.dir, filename)):
        fwords = npzfilename[:-len('.npz')].split('__')
        print npzfilename
        try: f = n.load(npzfilename)
        except(IOError): continue
        if fwords[1] == 'info':
            times = f['times']
            _afreqs = f['freqs']
            scores = f['scores']
            SHAPE = times.shape + _afreqs.shape
            filter = gen_filter(SHAPE, opts.dw, opts.drw)
            filter_take = n.where(filter)
            def from_coeffs(c):
                d = n.zeros(SHAPE, dtype=n.complex)
                d[filter_take] = c
                return fft2(d) / d.size
        else:
            k = fwords[1]
            srcs[k] = {}
            for i in f.files: srcs[k][int(i)] = f[i]
    best_score = scores.min()
    argclose = n.where(scores < best_score + 2*EXP_NOISE)[0]
    print len(argclose)
    print 'Using Score:', best_score
    srcant = {}
    for k in srcs:
        print k
        srcant[k] = {}
        for i in srcs[k]:
            _ant, _wgt = 0, 0
            for iter in argclose:
                w = n.exp((best_score - scores[iter]) / EXP_NOISE)
                _wgt += w
                _ant += from_coeffs(srcs[k][i][iter]) * w
            srcant[k][i] = _ant / _wgt
        if not srcdata.has_key(k): srcdata[k] = {}
        d = {}
        for i in srcant.get(k,{}):
          for j in srcant.get(k,{}):
            if j <= i: continue
            ai = srcant[k][i]
            aj = srcant[k][j]
            d[a.miriad.ij2bl(i,j)] = ai * n.conj(aj)
        flag = False
        for bl in d:
            srcdata[k][bl] = d[bl]
            flag = True
        
    srcs = srcs.keys()
    fq2fq = n.array([n.argmin(n.abs(_afreqs-f)).squeeze() for f in aa.get_afreqs()])
    #print fq2fq

    srclist = []
    for src in srcs:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)

    curtime, cnt = None, None
    def mfunc(uv, p, d, f):
        global curtime, cnt
        crd,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
            cat.compute(aa)
            cnt = n.argmin(n.abs(times-t))
            print t, cnt
        sd = 0
        for k in cat:
            try:
                try:
                    amp = srcdata[k][bl][cnt]
                    amp = amp.take(fq2fq)
                except(IndexError,KeyError): amp = 0
                phs = n.conj(aa.gen_phs(cat[k],i,j, srcshape=cat[k].srcshape, resolve_src=True))
            except(KeyError,IndexError,a.phs.PointingError): continue
            sd += amp * phs
            #try: sd += srcest_bl[k][bl][cnt].take(fq2fq) * phs / n.abs(phs).clip(1., n.Inf)
            #except(KeyError,IndexError): pass
        #try: sd += xtalk[bl].take(fq2fq)
        #except(KeyError): pass
        d -= sd * aa.passband(i,j)
        if i != j: f |= n.where(n.abs(d) > 1e0, 1, 0)
        return p, n.where(f, 0, d), f
        
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='RMNPZ: removed sources %s\n' % (str(srcs)))
