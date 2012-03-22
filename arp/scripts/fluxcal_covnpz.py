#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

EXP_NOISE = 1.

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

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
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
    for k in srcs:
        print k
        if not srcdata.has_key(k): srcdata[k] = {}
        flag = False
        for i in srcs[k]:
            _ant, _wgt = 0, 0
            for iter in argclose:
                w = n.exp((best_score - scores[iter]) / EXP_NOISE)
                _wgt += w
                _ant += from_coeffs(srcs[k][i][iter]) * w
            srcdata[k][i] = srcdata[k].get(i,[]) + [_ant / _wgt]
            flag = True
        if flag: srctimes[k] = srctimes.get(k,[]) + [times]
for k in srcdata:
    srctimes[k] = n.concatenate(srctimes[k], axis=0)
    for i in srcdata[k]:
        srcdata[k][i] = n.concatenate(srcdata[k][i], axis=0)
srcs = srcdata.keys(); srcs.sort()
if opts.cal != None:
    srclist = []
    for src in srcs:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)
    aa = a.cal.get_aa(opts.cal, afreqs)
else: cat = {}

#if 'cyg' in srcs: srcs = ['cyg'] + srcs
#norm=1
for cnt, k in enumerate(srcs):
  for i in srcdata[k]:
    #d,w = 0.,0.
    d = n.abs(srcdata[k][i])**2
    w = n.where(srcdata[k][i] == 0, 0, 1)
    t = srctimes[k]
    #d *= norm

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
    if True:
        bp = n.sqrt(n.where(spec == 0, 1, spec) / cat[k].jys)
        bp_poly = n.polyfit(afreqs, bp, deg=6)
        bp = n.polyval(bp_poly, afreqs)**2
    else: bp = aa.passband(i,i)
    spec /= bp
    #if cnt == 0 and k == 'cyg':
    #    norm.shape = (1,norm.size)
    #    continue
    valid = n.where(spec != 0, 1, 0)
    ind, flx = n.polyfit(n.log10(afreqs.compress(valid)/.150), n.log10(spec.compress(valid)), deg=1)
    
    print '%25s: ANT=%3d FLX=%6.1f IND=%+4.2f' % (k, i, 10**flx, n.round(ind,2))
    if True: print list(bp_poly)
    color = colors[cnt%len(colors)]
    p.loglog(1e3*afreqs, spec, color+',', label=k+' %d'%i)
    p.loglog(1e3*afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), color+'-', label=k)
    p.xlim(1e3*afreqs[0], 1e3*afreqs[-1])
    p.ylim(10,1e5)

#p.subplot(211)
#p.legend(loc='best')
p.grid()
p.show()

