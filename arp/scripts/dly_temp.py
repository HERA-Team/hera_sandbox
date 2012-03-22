#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, max=True, drng=True, dec=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])

MIN_CH,MAX_CH,SUBBAND = 60,804,50
aa.select_chans(n.arange(MIN_CH, MAX_CH))
freqs = aa.get_afreqs()
UV_RES = 1.5
LST_RES = 2*n.pi/36
pb_poly = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]
lstbins = n.arange(0,2*n.pi, LST_RES)

def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(LST_RES)
    return u,v, lst
def jy2T(freqs_in_GHz):
    lam = a.const.c / (freqs_in_GHz * 1e9)
    pb = n.polyval(pb_poly, freqs_in_GHz)
    return 1e-23 * lam**2 / (2 * a.const.k * pb)

data = {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    times = []
    src = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            aa.set_jultime(t)
            # need to bin in lst to accumulate coherently...
            if True: lst = lstbins[n.argmin(n.abs(lstbins - aa.sidereal_time()))]
            else: lst = aa.sidereal_time()
            src = a.phs.RadioFixedBody(lst, aa.lat)
            src.compute(aa)
            times.append(t)
            bm_wgt = aa[0].bm_response(src.get_crds('top'), pol='y')**2
            bm_wgt = bm_wgt.flatten()
        bl = a.miriad.ij2bl(i,j)
        u,v,w = aa.gen_uvw(i,j, src=src)
        u,v = u.flatten(), v.flatten()
        bin = uv2bin(n.average(u), n.average(v), lst)
        d = d[MIN_CH:MAX_CH]
        val = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        # mask out overly short baselines adversely affected by xtalk removal
        val *= n.where(n.abs(u) < 12.5, 0, 1)
        d = aa.phs2src(d, src, i, j)
        d /= aa.passband(i,j)
        d *= jy2T(freqs)
        d *= val
        if n.average(val) < .5: continue
        _d = n.fft.ifft(d)
        ker = n.fft.ifft(val)
        gain = n.abs(ker[0])
        if True:  # This can reduce src sidelobes
            _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            _d = info['res']
            d = n.fft.fft(_d) * val
        _d *= n.sqrt(float(_d.size)) / gain
        # Weight by square of primary beam response in the direction we're projecting for.
        #d *= bm_wgt
        #w = bm_wgt**2 * val
        if not data.has_key(bin): data[bin] = []
        if False: data[bin].append(d)
        else: data[bin].append(n.concatenate([_d[_d.size/2:],_d[:_d.size/2]]))

nplots = len(data)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)

pk_tot_sum, pk_tot_wgt = 0, 0
for cnt, bin in enumerate(data):
    d = n.array(data[bin])
    if False:
        if True: d = n.log10(n.abs(d).clip(1e-10,n.Inf))
        else: d = n.angle(d)
        dmax, dmin = None, None
        p.subplot(d1,d2,cnt + 1)
        if not opts.max is None: dmax = opts.max
        else: dmax = d.max()
        if not opts.drng is None: dmin = dmax - opts.drng
        else: dmin = d.min()
        p.imshow(d, aspect='auto', interpolation='nearest', vmax=dmax, vmin=dmin)
        p.colorbar()
    else:
        pk_sum, pk_wgt = 0, 0
        for i in range(d.shape[0]):
            for j in range(i+1, d.shape[0]):
                pk_sum += d[i]*n.conj(d[j])
                pk_wgt += 1
        if pk_wgt == 0: continue
        pk_tot_sum += pk_sum
        pk_tot_wgt += pk_wgt
        #p.semilogy(1e6*n.abs(pk_sum / pk_wgt))

p.semilogy(1e6*n.abs(pk_tot_sum / pk_tot_wgt))

p.show()
