import aipy as a, numpy as n

def sim_src_Jy(fq_GHz, amp_Jy_150=.1, index=-1., dly_ns=100., phs_offset=1.):
    spec = amp_Jy_150 * (fq_GHz / .150)**index
    phs = phs_offset * n.exp(-2j*n.pi * fq_GHz.astype(n.complex64) * dly_ns)
    return spec * phs

def sim_srcs(aa, pol, cnt_100Jy_per_sr_Jy=.1, cnt_index=-2., avg_index=-1., 
        std_index=.25, lo_cutoff_Jy=.3, hi_cutoff_Jy=100., bl_ns=100., 
        std_ionref_arcmin=1.):
    fq_GHz = aa.get_afreqs()
    bins = 10**n.arange(n.log10(hi_cutoff_Jy), n.log10(lo_cutoff_Jy), -.3)
    spec = n.zeros_like(fq_GHz, dtype=n.complex64)
    umag = bl_ns * fq_GHz
    for i, flx in enumerate(bins[:-1]):
        flx_interval = bins[i] - bins[i+1]
        flx += flx_interval/2
        sr = 2*n.pi
        cnt = int(n.around(cnt_100Jy_per_sr_Jy * flx_interval * sr * (flx / 100.)**cnt_index))
        # I'm ignoring beam weighting on sky
        for j in xrange(cnt):
            index = avg_index + n.random.normal(scale=std_index)
            while True:
                x = n.sin(n.random.uniform(-n.pi/2, n.pi/2))
                y = n.sin(n.random.uniform(-n.pi/2, n.pi/2))
                z = x**2 + y**2
                if z <= 1:
                    z = n.sqrt(1-z)
                    break
            bm_resp = aa[0].bm_response((x,y,z), pol=pol[0]) * \
                n.conj(aa[0].bm_response((x,y,z), pol=pol[1]))
            bm_resp = bm_resp.flatten()
            #bm_resp = n.exp(-(x**2+y**2) / (2*20*a.ephem.degree)**2)
            dly = x * bl_ns
            dth = n.random.normal(scale=std_ionref_arcmin) * a.ephem.arcminute
            dw = dth * umag * (fq_GHz / .150)**-2
            o = n.random.uniform()
            phs = n.exp(-2j*n.pi*(dw + o))
            spec += bm_resp * sim_src_Jy(fq_GHz, flx, index=index, dly_ns=dly, phs_offset=phs)
    return spec

def sim_sync(fq_GHz, hmap, aa, bl_ns=100., mfreq=.408, pol=None, jd=2455746.5):
    BL_SCALAR = 4.4
    SDF0 = fq_GHz[1] - fq_GHz[0]
    NCHAN = fq_GHz.size * 16
    MFREQ = n.average(fq_GHz)
    px = n.arange(hmap.npix())
    px_area = 4*n.pi / hmap.npix()
    if True:
        tx,ty,tz = hmap.px2crd(px)
        valid = n.where(tz > 0)
        tx,ty,tz = tx[valid], ty[valid], tz[valid]
    else:
        im = a.img.Img(size=400, res=.4)
        tx,ty,tz = im.get_top(center=(500,500))
        tx,ty,tz = tx.filled(0).flatten(), ty.filled(0).flatten(), tz.filled(-1).flatten()
    top = n.array([tx, ty, tz])

    if pol != None:
        afreqs = aa.get_afreqs()
        ch = n.argmin(n.abs(afreqs - MFREQ))
        aa.select_chans(n.array([ch]))
        p1,p2 = pol
        bm_resp = aa[0].bm_response((tx,ty,tz), pol=p1) * \
            n.conj(aa[0].bm_response((tx,ty,tz), pol=p2))
        bm_resp = bm_resp.flatten()
        aa.select_chans()
    else: bm_resp = 1

    bx,by,bz = 0.,0.,0.
    try: bx,by,bz = bl_ns
    except(ValueError):
        try: bx,by = bl_ns
        except(ValueError): bx = bl_ns
    BL_LEN = n.sqrt(bx**2 + by**2 + bz**2) * BL_SCALAR
    SDF = 2**n.floor(n.log2(1 / BL_LEN / SDF0)) * SDF0
    BW = NCHAN * SDF
    bins = n.fft.fftfreq(NCHAN, SDF)
    #bins = n.concatenate([bins[bins.size/2:], bins[:bins.size/2]])
    dbin = bins[1] - bins[0]
    tau = bx*tx + by*ty + bz*tz
    taubin = n.around(tau / dbin).astype(n.int)
    taures = tau - taubin * dbin
    phs_res = n.exp(-2j*n.pi*taures.astype(n.complex) * MFREQ)
    taubin.shape += (1,)

    aa.set_jultime(jd)

    # Precessing
    m_precess = a.coord.convert_m('eq','eq', oepoch=aa.epoch)
    m = n.linalg.inv(n.dot(aa.eq2top_m, m_precess)) # takes top to J2000
    ex,ey,ez = n.dot(m, top)
    d = hmap[ex,ey,ez] * bm_resp
    if False:
        import pylab as p
        d.shape = (1000,1000)
        p.imshow(n.abs(d))
        p.colorbar()
        p.show()

    # Misc preparations
    freq = n.fft.fftfreq(bins.size, dbin)
    freq = n.where(freq < 0, freq + BW, freq)
    while n.all(freq < MFREQ): freq += BW
    sdf = freq[1] - freq[0]
    window1 = n.where(freq + sdf/2 >= .05, 1., 0) * n.where(freq + sdf/2 < .25, 1., 0)
    freq_pad = n.arange(.05,.25, SDF0)
    sdf = freq_pad[1] - freq_pad[0]
    window2 = n.where(freq_pad + sdf/2 >= .1, 1., 0) * n.where(freq_pad + sdf/2 < .2, 1., 0)
    # Sub-delay bin phasing
    d_bl = d * phs_res

    # Binning
    hist = n.zeros(bins.size, dtype=n.complex)
    a.utils.add2array(hist, taubin, d_bl)

    # Computing coarse spectrum
    #hist = n.concatenate([hist[hist.size/2:], hist[:hist.size/2]])
    spec = n.fft.fft(hist)
    spec = spec.compress(window1)
    freq = freq.compress(window1)
    if False:
        import pylab as p
        p.subplot(121)
        p.plot(bins, hist)
        p.subplot(122)
        p.plot(freq, spec.real)
        #p.show()

    # Computing fine spectrum
    # window is necessary here to prevent high-freq band edges from coupling
    # into delay spectrum model, and then back out into the final spectrum.
    w = a.dsp.gen_window(spec.size, window='blackman-harris')
    dspec = n.fft.ifft(spec*w)
    dspec, info = a.deconv.clean(dspec, n.fft.ifft(w), tol=1e-9, stop_if_div=False)
    dspec_pad = n.zeros(2*fq_GHz.size, dtype=n.complex)
    dspec_pad[:dspec.size/2] = dspec[:dspec.size/2]
    dspec_pad[-dspec.size/2+1:] = dspec[-dspec.size/2+1:]
    spec_pad = n.fft.fft(dspec_pad)
    spec_pad = spec_pad.compress(window2)

    if False:
        import pylab as p
        p.subplot(121)
        p.plot(n.fft.fftfreq(freq.size, freq[1]-freq[0]), n.abs(dspec))
        p.subplot(122)
        p.plot(freq_pad.compress(window2), spec_pad.real)
        p.subplot(121)
        ddspec = n.fft.ifft(spec_pad)
        freq_pad = freq_pad.compress(window2)
        bins = n.fft.fftfreq(freq_pad.size, freq_pad[1]-freq_pad[0])
        p.plot(bins, n.abs(ddspec))
        p.show()

    # Apply synchrotron spectral index
    spec_pad *= (fq_GHz / mfreq)**-2.5
    #spec_pad *= (.150 / mfreq)**-2.5
    # Convert to Jy
    spec_pad *= 2 * a.const.k / (a.const.c/fq_GHz/1e9)**2 * px_area / 1e-23
    #spec_pad *= 2 * a.const.k / 200.**2 * px_area / 1e-23
    return spec_pad

