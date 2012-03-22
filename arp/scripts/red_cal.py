#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P

#filename = 'zen.2455600.70036.uvcmbRx'
filename = 'zen.2455600.74212.uvcmbRx'
uv = a.miriad.UV(filename)
fqs = n.arange(uv['nchan'])*uv['sdf'] + uv['sfreq']
del(uv)
#CH0,CH1 = 600,750
WIDTH = 150 # channels, ~7 MHz
CH0S = range(WIDTH,fqs.size-2*WIDTH,WIDTH/2)
#CH0,CH1 = 0,fqs.size
aa = a.cal.get_aa('pgb600_v001', fqs)
USE_VERTICAL = False
PLOT = True
VERBOSE = False
colors = 'kbgrycm'
Omega = 0.76
BL_LEN = 40e2
BL_HORIZON = 1.1
NOISE_WHITEN = True
#afq = n.average(fqs[CH0:CH1])
#UMAG = BL_LEN / (a.const.c / (1e9*afq))
#redshift = C.pspec.F21 / n.average(afq) - 1
#kpr = UMAG * C.pspec.dk_du(redshift)
#dk_deta = C.pspec.dk_deta(redshift)
#X2Y = C.pspec.X2Y(redshift)
#jy2T = C.pspec.jy2T(fqs[CH0:CH1])
#B = fqs[CH1] - fqs[CH0]
bl_sets = ['3_11' + ',3_19,4_20,5_21' + ',4_12,5_13,6_14,6_22',
           #'3_11,4_12,5_13,6_14,3_19,4_20,5_21,6_22,1_9,1_25,2_26,7_31',
           '3_12,4_13,5_14,4_19,5_20,6_21',
           '3_13,4_14,5_19,6_20',
           '3_14,6_19', # These have an add'l phase diff from one another
           '4_11,5_12,6_13,3_20,4_21,5_22',
           '5_11,6_12,3_21,4_22',
           '6_11,3_22',
           '11_19,12_20,13_21,14_22',
           '12_19,13_20,14_21',
           '13_19,14_20',
           '11_20,12_21,13_22',
           '11_21,12_22',
]
if USE_VERTICAL:
    bl_sets += [
           '3_4,4_5,5_6,11_12,12_13,13_14,19_20,20_21,21_22',
           '3_5,4_6,11_13,12_14,19_21,20_22',
           '3_6,11_14,19_22',
    ]
conjugate = '3_12,4_13,5_14,4_19,5_20,6_21'
conjugate += ',3_13,4_14,5_19,6_20'
conjugate += ',3_14,6_19'
conjugate = [a.miriad.ij2bl(int(word.split('_')[0]),int(word.split('_')[1])) for word in conjugate.split(',')]
bl_sets = [bl_sets[0]] # just do one set of redundant bls
t,d,f = C.arp.get_dict_of_uv_data([filename], ','.join(bl_sets), 'yy')
aa.set_jultime(.5*(t[0] + t[-1]))
print 'Time:', n.around((t[-1]-t[0]) / (3600*a.ephem.second),2), 'hrs @', aa.sidereal_time()
src = a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat)
src.compute(aa)
print 'Zenith:', src.ra, src.dec
bl_sets = [[a.miriad.ij2bl(int(word.split('_')[0]),int(word.split('_')[1])) 
        for word in bllist.split(',')]
            for bllist in bl_sets]
for bl in d:
    i,j = a.miriad.bl2ij(bl)
    if j in range(16,24) and not i in range(16,24): d[bl] = n.conjugate(d[bl])
    if bl in conjugate: d[bl] = n.conjugate(d[bl])

v = {}
for k in f: v[k] = n.logical_not(f[k])
for k in d: d[k] *= v[k]

# Calibrate all baselines to the first in the list & flatten fringes
pspec_sum,pspec_wgt = {},{}
pspec_fg_sum,pspec_fg_wgt = {},{}
for bltype,bls in enumerate(bl_sets):
    color = colors[bltype % len(colors)]
    pspec_sum[bltype] = {}; pspec_wgt[bltype] = {}
    pspec_fg_sum[bltype] = {}; pspec_fg_wgt[bltype] = {}
    cnt, bl0 = 0, bls[0]
    i0,j0 = a.miriad.bl2ij(bl0)
    if j0 in range(16,24) and not i0 in range(16,24): i0,j0 = j0,i0
    phs2zen = []
    for _t in t:
        aa.set_jultime(_t)
        src.compute(aa)
        phs2zen.append(aa.gen_phs(src, i0, j0))
        #x,y,z = aa.get_baseline(i0,j0, src=src)
        #u,v = x*.2, y*.2
        #print u,v#, aa.gen_uvw(i0,j0,src)[:2,0,-1]
    phs2zen = n.array(phs2zen)
    gain_factor = .01  # approximate Jy conv for now
    for bl in bls: d[bl] *= phs2zen * gain_factor
    #print n.average(n.abs(d[bl0]))
    window = a.dsp.gen_window(fqs.size, window='kaiser3')
    if False: # Fit delay of first baseline assuming symmetry around 0
        ds0 = d[bl0].sum(axis=0)
        vs0 = v[bl0].sum(axis=0)
        ds0 = ds0 / n.where(vs0 > 0, vs0, 1)
        tau,dtau = 0,0
        for j in range(4):
            ds0 *= n.exp(-1j*(fqs*dtau))
            tau += dtau
            _ds0 = n.abs(n.fft.ifft(window * ds0))
            noise_est = n.sqrt(n.average(n.abs(_ds0[_ds0.size/4:3*_ds0.size/4])**2))
            #print noise_est, _ds.max()
            _ds0 = n.where(_ds0 > 10*noise_est, _ds0, 0)
            _ds0 = _ds0.clip(0, 10*noise_est)
            taus = n.fft.fftfreq(ds0.size, d=fqs[1]-fqs[0])
            dtau = -2*n.pi*n.sum(_ds0 * taus) / n.sum(_ds0)
            #print dtau
        fixphs = n.exp(-1j*fqs*tau)
        fixphs.shape = (1,fixphs.size)
        d[bl0] *= fixphs
    
    # Remove gain/phase difference
    for bl1 in bls[cnt+1:]:
        i1,j1 = a.miriad.bl2ij(bl1)
        gain = (n.abs(d[bl1]).sum() / v[bl1].sum()) / (n.abs(d[bl0]).sum() / v[bl0].sum())
        d[bl1] /= gain
        # Permute indices to reflect conjugation above
        if j1 in range(16,24): i1,j1 = j1,i1
        # Compute measured values
        tau,off,dtau,doff = 0,0,0,0
        dij01 = d[bl0] * n.conj(d[bl1])
        dij01_sum = n.sum(dij01,axis=0)
        dij01_wgt = n.sum(v[bl0]*v[bl1],axis=0)
        for j in range(10):
            dij01_sum *= n.exp(-1j*(fqs*dtau+doff))
            tau += dtau; off += doff
            _phs = n.abs(n.fft.fft(dij01_sum))
            _wgt = n.abs(n.fft.fft(dij01_wgt))
            _phs,info = a.deconv.clean(_phs, _wgt, tol=1e-4)
            mx = n.argmax(_phs)
            if mx == 0:
                # Fine-tune calibration with linear fit
                valid = n.where(dij01_wgt > dij01_wgt.max()/2, 1, 0)
                fqs_val = fqs.compress(valid)
                dly = n.real(n.log(dij01_sum.compress(valid))/1j) # This doesn't weight data
                dtau,doff = n.polyfit(fqs_val, dly, deg=1)
            else:
                # Pull out an integral number of phase wraps
                if mx > _phs.size/2: mx -= _phs.size
                dtau,doff = 2*n.pi*mx / (fqs[-1] - fqs[0]), 0
        off %= 2*n.pi
        print ((i0,j0),(i1,j1)),':',n.around([tau,off], 3), n.around(gain,3)
        d[bl1] *= n.exp(1j*(fqs*tau+off)) # +1j insread of -1j b/c conj above
    
    # On with computing power spectrum
    dly,wgt = {}, {}
    dly_fg, wgt_fg = {}, {}
    bl_vs_t = {}
    for bl in bls:
        ds = d[bl].sum(axis=0)
        vs = v[bl].sum(axis=0)
        
        if True:
            gain = n.average(window*vs)
            _d = n.fft.ifft(window*ds) 
            _w = n.fft.ifft(window*vs)
        else:
            gain = n.average(vs)
            _d = n.fft.ifft(ds) 
            _w = n.fft.ifft(vs)
        c_d,info = a.deconv.clean(_d, _w, tol=1e-5, verbose=VERBOSE, stop_if_div=False, maxiter=200)
        taus = n.fft.fftfreq(c_d.size, d=fqs[1]-fqs[0])
        # Restrict clean components to within 10% of horizon length; P.grid()
        c_d = n.where(n.abs(taus) <= BL_HORIZON * BL_LEN / a.const.len_ns, c_d, 0)
        _c_d = n.fft.fft(c_d)
        ds -= vs * _c_d
        # Remove smooth model from original (unsummed) data to see how residuals beat down vs. time
        d[bl] -= v[bl] * _c_d.reshape(1,_c_d.size)
        if NOISE_WHITEN:
            whitener = n.random.normal(size=ds.size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=ds.size))
            #ds += 400 * (fqs/.150)**-2.5 * n.sqrt(vs.max()-vs) * whitener
            ds += 800 * (fqs/.150)**-2.5 * n.sqrt(vs.max()-vs) * whitener
            # Reweight vs to reflect that remainders have been filled in with noise
            vs = vs.max() * n.ones_like(vs)
            # Whiten the noise in the original (unsummed) data as well
            whitener = n.random.normal(size=d[bl].size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=d[bl].size))
            whitener.shape = d[bl].shape
            d[bl] += 800 * (fqs.reshape(1,fqs.size)/.150)**-2.5 * n.sqrt(v[bl].max()-v[bl]) * whitener
            v[bl] += v[bl].max() * n.ones_like(vs)
        if PLOT:
            P.subplot(241)
            P.plot(fqs, n.real(ds/n.where(vs>0,vs,1)), 'k,', alpha=0.2)
            P.plot(fqs, n.imag(ds/n.where(vs>0,vs,1)), 'r,', alpha=0.2)
            P.subplot(245)
            P.plot(fqs, n.angle(ds/n.where(vs>0,vs,1)), color+',', alpha=0.2)
            sums,wgts = [],[]
            for cnt in range(10):
                s,w = 0,0
                order = range(d[bl].shape[0])
                n.random.shuffle(order)
                for i in order:
                    s += d[bl][i]
                    #w += v[bl][i]
                    w += 1
                    sums.append(s.copy())
                    #wgts.append(w.copy())
                    wgts.append(w)
            print n.array(sums).shape, n.array(wgts).shape
            slope = n.polyfit(n.log(wgts), n.log(n.abs(sums)), deg=1)[0]
            P.subplot(242)
            P.plot(fqs, slope,',',alpha=0.5)
                
            
        for ch0 in CH0S:
            ch1 = ch0 + WIDTH
            if not dly.has_key(bl):
                dly[bl],wgt[bl] = {},{}
                dly_fg[bl],wgt_fg[bl] = {},{}
            if False:
                small_window = a.dsp.gen_window(WIDTH, window='kaiser3')
                gain = n.average(small_window*vs[ch0:ch1])
                dly[bl][ch0] = n.fft.ifft(small_window*ds[ch0:ch1])
                dly_fg[bl][ch0] = n.fft.ifft(small_window*_c_d[ch0:ch1])
            elif False:
                gain = n.average(vs[ch0:ch1])
                dly[bl][ch0] = n.fft.ifft(ds[ch0:ch1])
                dly_fg[bl][ch0] = n.fft.ifft(_c_d[ch0:ch1])
            else:
                gain = n.average(vs[ch0:ch1])
                dly[bl][ch0] = C.pfb.pfb(ds[ch0-WIDTH:ch1+WIDTH], window='kaiser3', taps=3, fft=n.fft.ifft)
                dly_fg[bl][ch0] = C.pfb.pfb(_c_d[ch0-WIDTH:ch1+WIDTH], window='kaiser3', taps=3, fft=n.fft.ifft)
            #dly[bl] = n.fft.ifft(ds-vs*n.fft.fft(c_d))
            wgt[bl][ch0] = gain
            wgt_fg[bl][ch0] = 1 
            #wgt[bl] = n.average(vs)
            if PLOT:
                P.subplot(241)
                P.plot([fqs[ch0],fqs[ch0]], [-1e6,1e6], 'g')
                P.plot([fqs[ch1],fqs[ch1]], [-1e6,1e6], 'g')
                #P.plot(fqs, n.abs(ds/n.where(vs>0,vs,1))/n.where(n.abs(ds0)>0, n.abs(ds0), 1), color+',', alpha=0.2)
                P.subplot(245)
                P.plot([fqs[ch0],fqs[ch0]], [-1e6,1e6], 'g')
                P.plot([fqs[ch1],fqs[ch1]], [-1e6,1e6], 'g')
    
    pspec_vs_bl_sum = {bltype:{}}
    pspec_vs_bl_wgt = {bltype:{}}
    # Sum up power spectrum
    for cnt,bl0 in enumerate(bls):
        for bl1 in bls[cnt+1:]:
            for ch0 in dly[bl0]:
                ch1 = ch0 + WIDTH
                jy2T = C.pspec.jy2T(fqs[ch0:ch1])
                if not pspec_sum[bltype].has_key(ch0):
                    pspec_sum[bltype][ch0] = 0; pspec_wgt[bltype][ch0] = 0
                    pspec_fg_sum[bltype][ch0] = 0; pspec_fg_wgt[bltype][ch0] = 0
                    pspec_vs_bl_sum[bltype][ch0] = []; pspec_vs_bl_wgt[bltype][ch0] = []
                dly2 = dly[bl0][ch0] * n.conj(dly[bl1][ch0]) * jy2T**2
                wgt2 = wgt[bl0][ch0] * wgt[bl1][ch0]
                dly_fg2 = dly_fg[bl0][ch0] * n.conj(dly_fg[bl1][ch0]) * jy2T**2
                wgt_fg2 = wgt_fg[bl0][ch0] * wgt_fg[bl1][ch0]
                pspec_sum[bltype][ch0] += dly2
                pspec_fg_sum[bltype][ch0] += dly_fg2
                pspec_wgt[bltype][ch0] += wgt2
                pspec_fg_wgt[bltype][ch0] += wgt_fg2
                pspec_vs_bl_sum[bltype][ch0].append(dly2)
                pspec_vs_bl_wgt[bltype][ch0].append(wgt2)

    if PLOT:
        P.subplot(246)
        P.title('Allen Variance Index vs. BL')
        for ch0 in pspec_vs_bl_sum[bltype]:
            sums,wgts = [],[]
            # compute several different summing orders to get sense of how noise beats down with adding
            for cnt in range(25):
                s,w = 0,0
                order = range(len(pspec_vs_bl_sum[bltype][ch0]))
                #print 'Working on %d baselines' % (order[-1]+1)
                n.random.shuffle(order)
                for i in order:
                    s += pspec_vs_bl_sum[bltype][ch0][i]
                    w += pspec_vs_bl_wgt[bltype][ch0][i]
                    sums.append(s.copy())
                    wgts.append(w.copy())
            slope = n.polyfit(n.log(wgts), n.log(n.abs(sums)), deg=1)[0]
            print n.array(sums).shape, n.array(wgts).shape
            P.plot(slope,',',alpha=0.5)
        
    afqs = []
    ks_all_ch, pspec_all_ch = [], []
    pspec_fg_all_ch = []
    cosmo_scalars = []
    for ch0 in CH0S:
        ch1 = ch0 + WIDTH
        afq = n.average(fqs[ch0:ch1])
        umag = BL_LEN / (a.const.c / (1e9*afq))
        redshift = C.pspec.F21 / n.average(afq) - 1
        kpr = umag * C.pspec.dk_du(redshift)
        dk_deta = C.pspec.dk_deta(redshift)
        x2y = C.pspec.X2Y(redshift)
        B = fqs[ch1] - fqs[ch0]
        print '-'*60
        print 'Channels %d-%d' % (ch0, ch0+WIDTH)
        print 'Umag:', umag, 'kpr:', kpr
        print 'Bandwidth:', n.around(B, 4), 'GHz @', afq
        print 'dk_deta:', dk_deta
        if PLOT:
            #P.subplot(234)
            plt_dly = n.sqrt(n.abs(pspec_sum[bltype][ch0] / pspec_wgt[bltype][ch0]))
            plt_dly_fg = n.sqrt(n.abs(pspec_fg_sum[bltype][ch0] / pspec_fg_wgt[bltype][ch0]))
            #plt_dly = n.concatenate([plt_dly[-100:], plt_dly[:101]])
            #P.semilogy(n.abs(n.arange(-100,101)), plt_dly, color+'-')
            taus = n.fft.fftfreq(plt_dly.size, d=fqs[1]-fqs[0])
            #P.semilogy(n.abs(taus), plt_dly, color+'-')
            #P.semilogy(n.abs(taus), plt_dly_fg, color+':')
            #P.semilogy([BL_HORIZON*BL_LEN/a.const.len_ns, BL_HORIZON*BL_LEN/a.const.len_ns], [1e-10,1e10], 'g')
            #P.xlim(0,taus[min(100,taus.size/2)])
            kpls = taus * dk_deta
            ks = n.sqrt(kpr**2 + kpls**2)
            print 'taus:', taus[:3]
            print 'kpls:', kpls[:3]
            print 'ks:', ks[:3]
            # The residual pspec
            pspec_fold_sum = pspec_sum[bltype][ch0][:ks.size/2]
            pspec_fold_sum[1:] += pspec_sum[bltype][ch0][:ks.size/2:-1]
            pspec_fold_sum[0] *= 2
            pspec_fold_wgt = 2*pspec_wgt[bltype][ch0]
            pspec_fold = pspec_fold_sum / pspec_fold_wgt
            kbins, pspec_fold_bin = C.pspec.rebin_log(ks[2:ks.size/2], pspec_fold[2:], nbins=10)
            kbins = n.concatenate([ks[:2], kbins])
            pspec_fold = n.concatenate([pspec_fold[:2], pspec_fold_bin])
            print 'kbins:', kbins[:3]
            pspec_fold = n.abs(pspec_fold)
            # The foreground pspec
            pspec_fold_fg_sum = pspec_fg_sum[bltype][ch0][:ks.size/2]
            pspec_fold_fg_sum[1:] += pspec_fg_sum[bltype][ch0][:ks.size/2:-1]
            pspec_fold_fg_sum[0] *= 2
            pspec_fold_fg_wgt = 2*pspec_fg_wgt[bltype][ch0]
            pspec_fg_fold = pspec_fold_fg_sum / pspec_fold_fg_wgt
            junk, pspec_fg_fold_bin = C.pspec.rebin_log(ks[2:ks.size/2], pspec_fg_fold[2:], nbins=10)
            pspec_fg_fold = n.concatenate([pspec_fg_fold[:2], pspec_fg_fold_bin])
            pspec_fg_fold = n.abs(pspec_fg_fold)
            #
            kpl_skycutoff = BL_HORIZON * BL_LEN / a.const.len_ns * dk_deta
            k_skycutoff = n.sqrt(kpr**2 + kpl_skycutoff**2)
            P.subplot(243)
            P.loglog([k_skycutoff, k_skycutoff], [1e-3,1e20], 'g')
            P.loglog(kbins, pspec_fold, color+'.-')
            P.loglog(kbins, pspec_fg_fold, color+'o-')
            # Set scaling to cosmological units
            # has B in numerator b/c ifft above divides vis by Nchan and Jy definition divides vis by ch bw
            # after squaring, vis^2 must be mult by 
            cosmo_scalar = x2y * kbins**3 / (2*n.pi**2) * Omega * B
            cosmo_scalars.append(cosmo_scalar)
            ks_all_ch.append(kbins)
            afqs.append(n.ones_like(kbins)*afq)
            pspec_all_ch.append(pspec_fold)
            pspec_fg_all_ch.append(pspec_fg_fold)
            P.subplot(247)
            P.loglog(kbins, pspec_fold*cosmo_scalar, color+'.-')
            P.loglog(kbins, pspec_fg_fold*cosmo_scalar, color+'o-')
    ks_all_ch = n.array(ks_all_ch)
    afqs = n.array(afqs)
    pspec_all_ch = n.array(pspec_all_ch)
    pspec_fg_all_ch = n.array(pspec_fg_all_ch)
    cosmo_scalars = n.array(cosmo_scalars)
    P.subplot(244)
    P.title('$P(k)$')
    # Note that there are sidelobes in pspec_fg_all that may compromise EoR in the future
    P.contourf(n.log10(ks_all_ch), afqs, n.log10(n.abs(pspec_all_ch+pspec_fg_all_ch)).clip(0.01,5.99), n.arange(0,6.1,.1))
    P.colorbar(shrink=0.5)
    P.ylim(afqs[0,0],afqs[-1,-1])
    P.xlim(n.log10(ks_all_ch[0,0]), n.log10(ks_all_ch[-1,-1]))
    P.subplot(248)
    P.title('$k^3P(k)$')
    # Note that there are sidelobes in pspec_fg_all that may compromise EoR in the future
    P.contourf(n.log10(ks_all_ch), afqs, n.log10(n.abs(pspec_all_ch+pspec_fg_all_ch)*cosmo_scalars).clip(5.01,9.99), n.arange(5,10.1,.1))
    P.colorbar(shrink=0.5)
    P.ylim(afqs[0,0],afqs[-1,-1])
    P.xlim(n.log10(ks_all_ch[0,0]), n.log10(ks_all_ch[-1,-1]))
    
#P.subplot(224)
#cosmo_scalar = X2Y * kbins**3 / (2*n.pi**2) * Omega * B
#P.loglog(kbins, pspec_fold*cosmo_scalar, 'k-.')

P.subplot(241); P.xlim(fqs[0], fqs[-1]); P.ylim(-1000,1000)
P.subplot(242); P.xlim(fqs[0], fqs[-1])
P.subplot(245); P.xlim(fqs[0], fqs[-1]); P.ylim(-n.pi, n.pi)
#P.subplot(234); P.ylim(1e-1,1e3); P.grid()
P.subplot(243); P.ylim(1e-1,1e8); P.grid()
P.subplot(247); P.ylim(1e1,1e11); P.grid()

#P.subplot(223)
#plt_dly = n.sqrt(n.abs(pspec_sum / pspec_wgt))
#plt_dly = n.concatenate([plt_dly[-100:], plt_dly[:101]])
#P.semilogy(n.arange(-100,101), plt_dly, 'k-')

'''

        P.subplot(133)
        w = a.dsp.gen_window(d[bl1].shape[-1], window='kaiser3')
        _sd = n.fft.ifft(w*d[bl1].sum(axis=0)) / n.sqrt(n.average(sw**2))
        _sd = n.concatenate([_sd[-100:], _sd[:101]])
        P.semilogy(n.arange(-100,101), n.abs(_sd), 'k-')
        #P.legend(loc='best')
    d = sum/n.where(wgt > 0, wgt, 1)

    P.subplot(134)
    P.imshow(n.log10(n.abs(d)), aspect='auto', vmax=6, vmin=4)
    P.colorbar(shrink=.5)

    P.subplot(136)
    _d = n.fft.ifft(sum*w.reshape(1,w.size))
    _w = n.fft.ifft(wgt*w.reshape(1,w.size))
    for i in range(_d.shape[0]):
        c_d,info = a.deconv.clean(_d[i],_w[i], tol=1e-4)
        _d[i] = c_d + info['res']
    _d = n.concatenate([_d[:,-20:], _d[:,:21]], axis=1)
    P.imshow(n.log10(n.abs(_d)), aspect='auto', vmax=5, vmin=3)
    P.colorbar(shrink=.5)

    P.subplot(235)
    P.imshow(n.angle(sum), aspect='auto')
      
    P.subplot(231)
    sw = wgt.sum(axis=0)
    sd = sum.sum(axis=0) / n.where(sw > 0, sw, 1)
    P.semilogy(fqs, n.abs(sd), 'k-,')
    P.subplot(232)
    P.plot(fqs, n.angle(sd), ',')

    P.subplot(233)
    _sd = n.fft.ifft(w*sum.sum(axis=0)) / n.sqrt(n.average(sw**2))
    _sw = n.fft.ifft(w*wgt.sum(axis=0)) / n.sqrt(n.average(sw**2))
    c_sd,info = a.deconv.clean(_sd, _sw, tol=1e-4)
    c_sd += info['res']# / n.sqrt(n.average((w*wgt.sum(axis=0))**2))
    _sd = n.concatenate([_sd[-100:], _sd[:101]])
    P.semilogy(n.arange(-100,101), n.abs(_sd), 'k-,')
    _sw = n.concatenate([_sw[-100:], _sw[:101]])
    P.semilogy(n.arange(-100,101), n.abs(_sw), 'r-,')
    c_sd = n.concatenate([c_sd[-100:], c_sd[:101]])
    P.semilogy(n.arange(-100,101), n.abs(c_sd), 'g-,')

    P.show()


#m = n.array(m)
#measured = n.array(measured)
#
#print ants
#n.set_printoptions(threshold=n.nan)
#print m
#print
#print measured
#print
#x,res,rank,s = n.linalg.lstsq(m, measured)
#print x
#print s

#P.subplot(133)
##P.plot(fqs, n.where(d01_avg_abs > 0, n.log(d01_avg/d01_avg_abs)/(1j*fqs), 0).real)
##d01 = C.arp.clean_transform(d[bl0] * n.conj(d[bl1]), n.logical_or(f[bl0],f[bl1]), axis=1)
##d01 = a.img.recenter(d01, (0,d01.shape[1]/2))
##d01 = d01[:,d01.shape[1]/2-25:d01.shape[1]/2+25]
##C.arp.waterfall(d01, drng=1.5)
##C.arp.waterfall(d01, mode='phs', mx=n.pi, drng=2*n.pi)
##P.title('cross')
##P.colorbar(shrink=.5)
#
#C.arp.waterfall(d[bl1]*n.exp(-1j*fqs*tau), mode='phs', mx=n.pi, drng=2*n.pi)
##C.arp.waterfall(d1, drng=1.5)
#P.title('6-14/c')

'''
P.show()
