#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

TRIM = True
#EXP_NOISE = 10.
EXP_NOISE = 1.
#EXP_NOISE = .1

try:
    import fftw3
    print 'Using FFTW FFT'
    _fftw3_dat, _fftw3_fwd, _fftw3_rev = None, None, None
    def fft2(d):
        global _fftw3_fwd, _fftw3_dat
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        elif _fftw3_dat[:].shape != d:
            _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_fwd()
        return _fftw3_dat
    def ifft2(d):
        global _fftw3_rev, _fftw3_dat
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='backward', flags=['measure'])
        elif _fftw3_dat[:].shape != d:
            _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='backward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_rev()
        return _fftw3_dat
except(ImportError):
    print 'Using numpy FFT'
    fft2, ifft2 = n.fft.fft2, n.fft.ifft2

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True,src=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-q', '--quality', dest='quality', default=0., help='Cutoff for plotting a source.')
opts, args = o.parse_args(sys.argv[1:])

p.rcParams['legend.fontsize'] = 6

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

filegroups = {}
for cnt, filename in enumerate(args):
    basefile = filename.split('__')[0]
    filegroups[basefile] = filegroups.get(basefile, []) + [filename]
srcdata, srctimes = {}, {}
basefiles = filegroups.keys(); basefiles.sort()
for basefile in basefiles:
    filenames = filegroups[basefile]; filenames.sort(); filenames.reverse()
    srcs = {}
    for filename in filenames:
        fwords = filename[:-len('.npz')].split('__')
        print filename
        try: f = n.load(filename)
        except:
            print '    Load file failed'
            continue
        if fwords[1] == 'info':
            times = f['times']
            afreqs = f['freqs']
            scores = f['scores']
            SHAPE = times.shape + afreqs.shape
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
                _ant += srcs[k][i][iter] * w
            srcant[k][i] = from_coeffs(_ant / _wgt)
            if TRIM:
                trim = len(srcant[k][i]) / 4
                srcant[k][i] = srcant[k][i][trim:-trim]
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
            srcdata[k][bl] = srcdata[k].get(bl,[]) + [d[bl]]
            flag = True
        if flag:
            if TRIM: srctimes[k] = srctimes.get(k,[]) + [times[trim:-trim]]
            else: srctimes[k] = srctimes.get(k,[]) + [times]
for k in srcdata:
    srctimes[k] = n.concatenate(srctimes[k], axis=0)
    for bl in srcdata[k]:
        srcdata[k][bl] = n.concatenate(srcdata[k][bl], axis=0)
srcs = srcdata.keys(); srcs.sort()
if opts.cal != None:
    srclist = []
    for src in srcs:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    print srclist,[opts.cat]
    cat = a.cal.get_catalog(opts.cal, srcs=srclist,catalogs=[opts.cat])
    print cat.keys()
    aa = a.cal.get_aa(opts.cal, afreqs)
else: cat = {}

if 'cyg' in srcs: srcs = ['cyg'] + srcs
norm=1
for cnt, k in enumerate(srcs):
    d,w = 0.,0.
    for bl in srcdata[k]:
        d += srcdata[k][bl]
        w += n.where(srcdata[k][bl] == 0, 0, 1)
    #d /= w.clip(1,n.Inf)
    t = srctimes[k]
    #order = n.argsort(t)
    #d,t = d.take(order, axis=0), t.take(order)
    #I = 1
    #shape = (int(t.shape[0]/I), I)
    #ints = shape[0] * shape[1]
    #d,t = d[:ints], t[:ints]
    #d.shape,t.shape = shape + d.shape[1:], shape
    #d,t = n.average(d, axis=1), n.average(t, axis=1)
    d *= norm

    # Calculate beam response
    bm = []
    lsts = []
    for jd in t:
        aa.set_jultime(jd)
        lsts.append(aa.sidereal_time())
        cat[k].compute(aa)
        bm.append(aa[0].bm_response(cat[k].get_crds('top'), pol=opts.pol)**2)
    bm = n.array(bm).squeeze()
    #spec = n.sum(d*bm, axis=0)/n.sum(bm**2, axis=0)
    d_bm = n.sum(d*bm, axis=0)
    w_bm = n.sum(w*bm**2, axis=0)
    spec = d_bm / n.where(w_bm == 0, 1, w_bm)
    dsum = n.sum(n.abs(d), axis=1)
    wsum = n.sum(w, axis=1)
    vs_time = dsum / n.where(wsum == 0, 1, wsum)
    #bsum = n.sum(n.abs(cat[k].get_jys()*w*bm), axis=1)
    #bsum = n.sum(n.abs(w*bm), axis=1)
    bsum = n.sum(n.abs(100*(afreqs/.150)**-1*w*bm), axis=1)
    bm_vs_time = bsum / n.where(wsum == 0, 1, wsum)
    if cnt == 0 and k == 'cyg':
        norm = cat['cyg'].jys / n.where(spec == 0, 1, spec)
        norm.shape = (1,norm.size)
        continue
    valid = n.where(spec != 0, 1, 0)
    ind, flx = n.polyfit(n.log10(afreqs.compress(valid)/.150), n.log10(spec.compress(valid)), deg=1)
    
    d, w, bm = d.flatten(), w.flatten(), bm.flatten()
    valid = n.where(w == 0, 0, 1)
    d, bm = d.compress(valid), bm.compress(valid)
    q = n.average((d-n.average(d))*(bm - n.average(bm))) / n.std(d) / n.std(bm)
    #q = n.median((d-n.median(d))*(bm - n.median(bm))) / n.std(d) / n.std(bm)
    print '%25s: FLX=%6.1f IND=%+4.2f Q=%+4.2f' % (k, 10**flx, n.round(ind,2), n.round(q,2))
    #d /= bm
    #_f = flx[1:-1] - .5 * (flx[2:] + flx[:-2])
    #q = n.sum(n.abs(_f)) / n.sum(n.abs(flx[1:-1]))
    if q < opts.quality: continue
    color = colors[cnt%len(colors)]
    p.subplot(211)
    d.shape = (t.size,d.size/t.size)
    p.semilogy(t, n.average(n.abs(d), axis=1), color+',', label=k)
#    p.semilogy(t, n.abs(d), color+',', label=k)
    #p.semilogy(lsts, vs_time, color+',', label=k)
    #p.semilogy(lsts, bm_vs_time, color+'-', label=k)
    p.ylim(.1,1e5)

    p.subplot(212)
    p.loglog(afreqs, spec, color+',', label=k)
    #p.loglog(afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), color+'-', label=k)
    p.xlim(afreqs[0], afreqs[-1])
    p.ylim(10,1e5)

#p.subplot(211)
p.legend(loc='best')
p.show()

