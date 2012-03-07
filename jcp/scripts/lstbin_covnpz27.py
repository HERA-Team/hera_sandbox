#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a, math
import optparse
import qPickle

TRIM = True
#EXP_NOISE = 10.
EXP_NOISE = 1.
#EXP_NOISE = .1

try:
    import fftw3
    print 'Using FFTW FFT'
    _fftw3_dat, _fftw3_fwd, _fftw3_rev = None, None, None
    def fft2(d):
        global _fftw3_fwd, _fftw3_rev, _fftw3_dat
        if not _fftw3_dat is None and _fftw3_dat.shape != d.shape:
            _fftw3_fwd, _fftw3_rev, _fftw3_dat = None, None, None
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_fwd()
        return _fftw3_dat
    def ifft2(d):
        global _fftw3_fwd, _fftw3_rev, _fftw3_dat
        if not _fftw3_dat is None and _fftw3_dat.shape != d.shape:
            _fftw3_fwd, _fftw3_rev, _fftw3_dat = None, None, None
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
o.add_option('-q', '--quality', dest='quality', default=0., help='Cutoff for plotting a source.')
o.add_option('-a', '--ant', dest='ant', type=int,
    help='Which antenna to use in the plot. Default = All.')
o.add_option('-f', '--freq', dest='freq', type=float, default=.150,
    help='Frequency at which to evaluate the beam in GHz.  Default=.150')
o.add_option('--b', dest='balun', action='store_true',
    help='Plot the balun data.')
o.add_option('--c', dest='cable', action='store_true',
    help='Plot the cable data.')
opts, args = o.parse_args(sys.argv[1:])

#p.rcParams['legend.fontsize'] = 6

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
srcdata, srctimes, t_balun, t_cable = {}, {}, {}, {}
basefiles = filegroups.keys(); basefiles.sort()
antvis = {}
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
            cable = f['t_cable']
            balun = f['t_balun']
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
        if not antvis.has_key(k): antvis[k] = {}
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
                trim = len(srcant[k][i]) / 3
                srcant[k][i] = srcant[k][i][trim:-trim]
            if not antvis[k].has_key(i):
                antvis[k][i] = srcant[k][i]
            else:
                antvis[k][i] = n.append(antvis[k][i],srcant[k][i],axis=0)
                #print antvis[k][i].shape
        if TRIM:
            srctimes[k] = srctimes.get(k,[]) + [times[trim:-trim]]
            t_cable[k] = t_cable.get(k, []) + [cable[trim:-trim]]
            t_balun[k] = t_balun.get(k, []) + [balun[trim:-trim]]
        else:
            srctimes[k] = srctimes.get(k,[]) + [times]
            t_cable[k] = t_cable.get(k, []) + [cable]
            t_balun[k] = t_balun.get(k, []) + [balun]
for k in antvis:
    srctimes[k] = n.concatenate(srctimes[k], axis=0)
    t_cable[k] = n.concatenate(t_cable[k], axis=0)
    t_balun[k] = n.concatenate(t_balun[k], axis=0)
srcs = srcant.keys(); srcs.sort()
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

#if 'cyg' not in srcs: srcs = ['cyg'] + srcs
norm=1

vis_list = {}
cable_list = {}
cable_fit = {}
balun_list = {}
balun_fit = {}
ants = {}
for k in srcs:
    ants[k] = []
    vis_list[k] = {}
    cable_list[k] = {}
    cable_fit[k] = []
    balun_list[k] = {}
    balun_fit[k] = []
    for i in srcant[k]:
        if not vis_list[k].has_key(i):
            vis_list[k][i] = []
            cable_list[k][i] = []
            balun_list[k][i] = []
            ants[k].append(i)
        if i in [5,14]: continue
        for cnt, integ in enumerate(antvis[k][i]):
            #if cnt !=0:
            #    if t_balun[k][cnt] == _t_balun: continue
            #_t_balun = t_balun[k][cnt]
            aa.set_jultime(srctimes[k][cnt])
            cat.compute(aa)
            bm = aa[i].bm_response(cat[k].get_crds('top'),pol=opts.pol)
            gain = bm.squeeze()**2
            visint = ((integ*n.conj(integ))/(cat[k].jys*gain))
            vis_list[k][i].append(n.median(visint))
            cable_list[k][i].append(t_cable[k][cnt])
            #cable_fit[k].append(t_cable[k][cnt])
            balun_list[k][i].append(t_balun[k][cnt])
            #balun_fit[k].append(t_balun[k][cnt])

#t = n.histogram(t_cable['cas'],bins=2500)[0]
#p.plot(t, '.')

vis_fit = {}
residuals = {}
print srcs
for k in srcs:
    if opts.ant:
        print 'user defined antenna'
        ants[k] = [opts.ant]
    vis_fit[k] = []
    for i in ants[k]:
        vis_list[k][i] = vis_list[k][i]/n.median(vis_list[k][i])
        if len(vis_fit[k]) == 0:
            vis_fit[k] = vis_list[k][i]
            cable_fit[k] = cable_list[k][i]
            balun_fit[k] = balun_list[k][i]
        else:
            vis_fit[k] = n.append(vis_fit[k], vis_list[k][i])
            cable_fit[k] = n.append(cable_fit[k], cable_list[k][i])
            balun_fit[k] = n.append(balun_fit[k], balun_list[k][i])

    if opts.cable:
        m, b = n.polyfit(cable_fit[k],10*n.log10(vis_fit[k]),1)
        m, b = n.real(m), n.real(b)
        print m, b
        xplot = n.arange(n.floor(min(cable_fit[k])),n.ceil(max(cable_fit[k]))+1)
        p.title("Temperature Dependent Gains")
        p.subplot(121)
        for i in ants[k]:
            p.plot(cable_list[k][i], vis_list[k][i],'.')
        p.plot(xplot, (10**(((m*xplot) + b)/10)), 'black', linewidth=4)
        p.ylabel("Gain")
        p.xlabel("Cable Temperature (K)")
        #print len(vis_fit[k]), len(m*cable_fit[k])
        residuals[k] = n.real(sum((vis_fit[k]-(10**((m*cable_fit[k] + b)/10)))**2)/(len(vis_fit[k])))

    if opts.balun:
        m, b = n.polyfit(balun_fit[k],10*n.log10(vis_fit[k]),1)
        m, b = n.real(m), n.real(b)
        print m, b
        xplot = n.arange(n.floor(min(balun_fit[k])),n.ceil(max(balun_fit[k]))+1)
        p.subplot(122)
        for i in ants[k]:
             p.plot(balun_list[k][i], vis_list[k][i],'.')
        p.plot(xplot, (10**(((m*xplot) + b)/10)), 'black', linewidth=4)
        p.xlabel("PDA Temperature (K)")
        residuals[k] = n.real(sum((vis_fit[k]-(10**((m*balun_fit[k] + b)/10)))**2)/(len(vis_fit[k])))

    vis_fit[k] = vis_fit[k] - n.mean(vis_fit[k]) + 1.
    #print 'residuals', residuals[k]

    offset = n.ones_like(balun_fit[k])
    #X = n.matrix(n.column_stack((balun_fit[k],offset)))
    X = n.matrix(n.column_stack((balun_fit[k],cable_fit[k],offset)))
    Y = n.matrix((10 * n.log10(vis_fit[k]))).T.real
    a = (X.T * X)
    a = a.I
    b = X.T * Y
    H = n.dot(a,b)
    print "balun, cable, offset"
    print H
    H = H.tolist()

#residuals
#fit
    fit = H[0]*balun_fit[k] + H[1]*cable_fit[k] + H[2]
    #fit = fit - n.mean(fit)
    fit_residuals = n.real(sum((vis_fit[k] - (10**(fit/10)))**2)/(len(vis_fit[k])-2))
    print len(fit)

#polynomial model
    poly_model_b = 30.3573 - 0.02485*(balun_fit[k] - 273.) + 0.00010256*(balun_fit[k] - 273.)**2 - 0.000001979*(balun_fit[k] - 273.)**3
    poly_model_c = 15*(-0.72167 - 0.0032929*(cable_fit[k] - 273.) + 0.000078251*(cable_fit[k] - 273.)**2 - 0.0000013723*(cable_fit[k] - 273.)**3 + 0.0000000098601*(cable_fit[k] - 273.)**4)
    poly_model = poly_model_b + poly_model_c
    poly_model = poly_model - n.mean(poly_model)
    poly_model_residuals = n.real(sum((vis_fit[k] - (10**(poly_model/10)))**2)/(len(vis_fit[k])))

#linear model
    lin_model = n.real(-0.024*balun_fit[k] + -0.018*cable_fit[k])
    lin_model = lin_model - n.mean(lin_model)
    lin_model_residuals = n.real(sum((vis_fit[k] - (10**(lin_model/10)))**2)/(len(vis_fit[k])))

#null hypothesis
    null_model = n.real(0*balun_fit[k] + 0*cable_fit[k])
    null_model_residuals = n.real(sum((vis_fit[k] - (10**(null_model/10)))**2)/(len(vis_fit[k])))


    print 'fit'
    print fit_residuals
    print 'lin model'
    print lin_model_residuals
    print 'poly model'
    print poly_model_residuals
    print 'null hypothesis'
    print null_model_residuals

#p.plot(fit,lin_model,'.')
#p.plot(fit,poly_model,'.')
x = n.arange(-.4,.4,.001)
#p.plot(x,x)
#p.legend(("Linear", "Polynomial"),loc='upper left')
#p.xlabel("Gain Fluctuations from Fit")
#p.ylabel("Gain Fluctutions from Model")
#p.title("Model Comparisons")
p.show()
