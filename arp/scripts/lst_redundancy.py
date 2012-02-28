#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

LST_FINEBIN = 2e-3
LST_COARSEBIN = 1e-2
#LST_COARSEBIN = 1e-3
uv = a.miriad.UV(args[0])
NCHAN = uv['nchan']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
fqs = aa.get_afreqs()
del(uv)

dat,flg = {}, {}
for filename in args:
    print 'Reading', filename
    _dat,_flg = {}, {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Gather data from file
    curtime = None
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
        if n.all(f): continue
        #x,y,z = aa.get_baseline(i,j)
        x,y = 10*i,10*j
        bin = C.pspec.uv2bin(x,y,aa.sidereal_time(), lst_res=LST_FINEBIN)
        _dat[bin] = _dat.get(bin,0) + d
        _flg[bin] = _flg.get(bin,0) + n.logical_not(f).astype(n.int)
    for bin in _dat.keys():
        valid = _flg[bin]
        dat[bin] = dat.get(bin,[]) + [_dat[bin]/valid.clip(1,n.Inf)]
        flg[bin] = flg.get(bin,[]) + [valid.clip(0,1)]
        del(_dat[bin]); del(_flg[bin])

print 'Regridding to coarse lsts bins'
d00,d01,d11 = {},{},{}
w = {}
CH_STEP = 32
for bin in dat.keys():
    if len(dat[bin]) == 1:
        del(dat[bin])
        continue
    d0,d1 = dat[bin][0], dat[bin][-1]
    f0,f1 = flg[bin][0], flg[bin][-1]
    valid = f0 * f1
    d0 *= valid; d1 *= valid
    x,y,lst = C.pspec.bin2uv(bin, lst_res=LST_FINEBIN)
    bin = C.pspec.uv2bin(x,y,lst,lst_res=LST_COARSEBIN)
    d00[bin] = d00.get(bin,0) + n.abs(d0)**2
    d11[bin] = d11.get(bin,0) + n.abs(d1)**2
    d01[bin] = d01.get(bin,0) + d0 * n.conj(d1)
    w[bin] = w.get(bin,0) + valid

ants, lsts = {}, {}
for bin in d00:
    x,y,lst = C.pspec.bin2uv(bin,lst_res=LST_COARSEBIN)
    i,j = int(n.around(x/10)),int(n.around(y/10))
    ants[i] = ants[j] = None
    lsts[lst] = None
ants = ants.keys(); ants.sort()
lsts = lsts.keys(); lsts.sort()
bldat = {}
for cnt1,i in enumerate(ants):
  for j in ants[cnt1+1:]:
    bl = a.miriad.ij2bl(i,j)
    d = []
    for lst in lsts:
        bin = C.pspec.uv2bin(i*10,j*10,lst,lst_res=LST_COARSEBIN)
        _w = w[bin].sum()
        d.append([d00[bin].sum() / _w, d11[bin].sum() / _w, n.abs(d01[bin].sum()) / _w])
    bldat[bl] = n.array(d)
bldat['lst'] = n.array(lsts)
print 'Saving to', args[0]+'.lst.npz'
n.savez(args[0]+'.lst.npz', **bldat)
sys.exit(0)

print 'Self-calibrating gains'
g,g_sub,tau = {}, {}, {}
for cnt1,i in enumerate(ants):
  print cnt1+1, '/', len(ants)
  if not g.has_key(i): g[i] = {}
  for cnt2, j in enumerate(ants[cnt1+1:]):
    if not g.has_key(j): g[j] = {}
    cnt2 += cnt1 + 1
    for k in ants[cnt2+1:]:
      if not g.has_key(k): g[k] = {}
      for lst in lsts:
        bij = C.pspec.uv2bin(10*i,10*j,lst, lst_res=LST_COARSEBIN)
        bik = C.pspec.uv2bin(10*i,10*k,lst, lst_res=LST_COARSEBIN)
        bjk = C.pspec.uv2bin(10*j,10*k,lst, lst_res=LST_COARSEBIN)
        dij = n.sqrt(d00[bij] / d11[bij])
        dik = n.sqrt(d00[bik] / d11[bik])
        djk = n.sqrt(d00[bjk] / d11[bjk])
        #print i,j,k, dij, dik, djk
        g[i][lst] = g[i].get(lst,[]) + [dij * dik / djk]
        g[j][lst] = g[j].get(lst,[]) + [dij * djk / dik]
        g[k][lst] = g[k].get(lst,[]) + [dik * djk / dij]
avg = 0
_g = {}
for i in g:
    for lst in lsts: g[i][lst] = n.sqrt(n.average(g[i][lst]))
    _g[str(i)] = n.array([g[i][lst] for lst in lsts])
    avg += _g[str(i)]
avg /= len(g)
g = _g
g['lst'] = n.array(lsts)
print 'Saving to', args[0]+'.lst_redgain.npz'
n.savez(args[0]+'.lst_redgain.npz', **g)
sys.exit(0)

for ch in range(CH_STEP,NCHAN-CH_STEP,CH_STEP): g_sub[ch] = {}
for bin in d00:
    g[lst] = n.sqrt(d00[lst].sum() / d11[lst].sum())
    for ch in g_sub:
        _ch = ch+CH_STEP
        g_sub[ch][lst] = n.sqrt(d00[lst][ch:_ch].sum() / d11[lst][ch:_ch].sum())
    # Compute delay difference
    _tau,off,dtau,doff = 0,0,0,0
    dij01_sum = d01[lst]
    dij01_wgt = w[lst]
    for j in range(10):
        dij01_sum *= n.exp(-1j*(fqs*dtau+doff))
        _tau += dtau; off += doff
        _phs = n.abs(n.fft.fft(dij01_sum))
        _wgt = n.abs(n.fft.fft(dij01_wgt))
        _phs,info = a.deconv.clean(_phs, _wgt, tol=1e-3)
        mx = n.argmax(_phs)
        if mx == 0:
            # Fine-tune calibration with linear fit
            valid = n.where(dij01_wgt > dij01_wgt.max()/2., 1, 0)
            fqs_val = fqs.compress(valid)
            dly = n.real(n.log(dij01_sum.compress(valid))/1j) # This doesn't weight data
            dtau,doff = n.polyfit(fqs_val, dly, deg=1)
        else:
            # Pull out an integral number of phase wraps
            if mx > _phs.size/2: mx -= _phs.size
            dtau,doff = 2*n.pi*mx / (fqs[-1] - fqs[0]), 0
    off %= 2*n.pi
    tau[lst] = _tau
    print lst,':', g[lst], tau[lst]

lsts = g.keys(); lsts.sort()
g = n.array([g[lst] for lst in lsts])
tau = n.array([tau[lst] for lst in lsts])
for ch in g_sub: g_sub[ch] = n.array([g_sub[ch][lst] for lst in lsts])

p.subplot(211); p.plot(lsts, g, '-')
p.subplot(212); p.plot(lsts, tau, '-')
for ch in g_sub:
    continue
    label = str(int(n.around(fqs[ch+CH_STEP/2] * 1e3)))
    p.subplot(211); p.plot(lsts, g_sub[ch]/g, '-', label=label)
p.subplot(211)
p.ylabel('Gain [Day1/Day2]')
p.legend()
p.ylim(.5,1.5)
p.subplot(212)
p.ylabel('$\Delta$ Delay [ns]')
p.xlabel('LST')
p.ylim(-50,50)
p.show()
