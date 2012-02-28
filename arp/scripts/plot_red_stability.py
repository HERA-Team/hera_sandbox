#! /usr/bin/env python
'''
Does data acquisition for max-redundancy data, then calibrates phase based on phase-offset/delay model.
'''
import aipy as a, numpy as n, capo as C, pylab as P
import sys

DECONV = False

filenames = sys.argv[1:]
uv = a.miriad.UV(filenames[0])
fqs = n.arange(uv['nchan'])*uv['sdf'] + uv['sfreq']
del(uv)

t,d,f = C.arp.get_dict_of_uv_data(filenames, '3_11,4_12', 'yy', verbose=True)

conjugate = '3_12,4_13,5_14,4_19,5_20,6_21'
conjugate += ',3_13,4_14,5_19,6_20'
conjugate += ',3_14,6_19'
conjugate = [a.miriad.ij2bl(int(word.split('_')[0]),int(word.split('_')[1])) for word in conjugate.split(',')]
for bl in d:
    i,j = a.miriad.bl2ij(bl)
    if j in range(16,24) and not i in range(16,24): d[bl] = n.conjugate(d[bl])
    if bl in conjugate: d[bl] = n.conjugate(d[bl])

v = {}
for k in f: v[k] = n.logical_not(f[k])
for k in d: d[k] *= v[k]

window = a.dsp.gen_window(fqs.size, window='kaiser3')
# Calibrate all baselines to the first in the list & flatten fringes
bls = d.keys()
for cnt,bl0 in enumerate(bls):
    for bl1 in bls[cnt+1:]:
        # Compute measured values
        d01 = d[bl0] * n.conj(d[bl1]) * window
        v01 = v[bl0] * n.conj(v[bl1])
        sh = d[bl0].shape
        d0 = (d[bl0]*v01)[:sh[0]-(sh[0]%8),:sh[1]-(sh[1]%8)].reshape(sh[0]/8,8,sh[1]/8,8)
        d1 = (d[bl1]*v01)[:sh[0]-(sh[0]%8),:sh[1]-(sh[1]%8)].reshape(sh[0]/8,8,sh[1]/8,8)
        gain = n.abs(d0.sum(axis=3).sum(axis=1) / d1.sum(axis=3).sum(axis=1))
        v01 *= window
        val,dly,off,phs,times = [],[],[],[],[]
        NINT = 4
        for i in range(d01.shape[0]/NINT): # step through integrations
            #d0i,v0i = d[bl0][i], v[bl0][i]
            #d1i,v1i = d[bl1][i], v[bl1][i]
            d0i,v0i = d[bl0][NINT*i:NINT*(i+1)].sum(axis=0), v[bl0][NINT*i:NINT*(i+1)].sum(axis=0)
            d1i,v1i = d[bl1][NINT*i:NINT*(i+1)].sum(axis=0), v[bl1][NINT*i:NINT*(i+1)].sum(axis=0)
            d01i = d0i * n.conj(d1i)
            v01i = v0i * v1i
            if v01i.sum() < v01i.size*NINT/2:
                print i, v01i.sum(), '*'
                continue
            else: print i, v01i.sum(),
            times.append(t[i])
            val.append(n.logical_and(v0i,v1i))
            dlyi,offi,dtau,doff = 0,0,0,0
            if False:
                dlyi,offi = 87,0
                d01i *= n.exp(-1j*(fqs*dlyi+offi))
            else:
                mx = 1
                for j in range(10):  # number of steps to assure convergence?
                    d01i *= n.exp(-1j*(fqs*dtau+doff))
                    dlyi += dtau; offi += doff
                    if mx != 0: # If mx is already 0, then it should be so in perpetuity
                        _dat = n.abs(n.fft.fft(d01i))
                        _wgt = n.abs(n.fft.fft(v01i))
                        if DECONV: _dat,info = a.deconv.clean(_dat, _wgt, tol=1e-3)
                        mx = n.argmax(_dat)
                    if mx == 0:
                        # Fine-tune calibration with linear fit
                        fqs_val = fqs.compress(v01i > 0)
                        dat = n.real(n.log(d01i.compress(v01i > 0))/1j) # This doesn't weight data
                        #dtau,doff = n.polyfit(fqs_val, dat, deg=1)
                        dtau,doff = n.median(dat/fqs_val),0
                    else:
                        # Pull out an integral number of phase wraps
                        if mx > _dat.size/2: mx -= _dat.size
                        dtau,doff = 2*n.pi*mx / (fqs[-1] - fqs[0]), 0
            #dat = n.real(n.log(d01i)/1j) # This doesn't weight data
            #dat = n.angle(d01i) # This doesn't weight data
            dat = n.angle(d0i*n.conj(d1i)) # This doesn't weight data
            P.plot(fqs, dat, 'k,', alpha=0.05)
            P.plot(fqs, ((fqs*dlyi+offi+n.pi)%(2*n.pi))-n.pi, 'g-', alpha=0.05)
            print dlyi, offi
            phs.append(d01i)
            dly.append(dlyi)
            off.append(offi)
        P.show()
        phs = n.array(phs)
        phs = phs[:phs.shape[0]-(phs.shape[0]%8),:phs.shape[1]-(phs.shape[1]%8)]
        val = n.array(val)
        val = val[:val.shape[0]-(val.shape[0]%8),:val.shape[1]-(val.shape[1]%8)]
        phs.shape = (phs.shape[0]/8, 8, phs.shape[1]/8, 8)
        val.shape = (val.shape[0]/8, 8, val.shape[1]/8, 8)
        phs = n.angle(phs.sum(axis=3).sum(axis=1) / val.sum(axis=3).sum(axis=1))
        dly = n.array(dly)
        off = ((n.array(off) + n.pi) % (2*n.pi)) - n.pi
        times = n.array(times)
        #print ((i0,j0),(i1,j1)),':',n.around([tau,off], 3), n.around(gain,3)
        #P.subplot(311); P.plot(times,gain,',',alpha=0.2)
        P.subplot(221); P.imshow(n.log10(gain), vmax=.3, vmin=-.3, aspect='auto', interpolation='nearest'); P.colorbar(shrink=.5)
        P.subplot(222); P.imshow(phs, vmax=1, vmin=-1, aspect='auto', interpolation='nearest'); P.colorbar(shrink=.5)
        P.subplot(223); P.plot(times,dly,',')
        P.subplot(224); P.plot(times,off,',')
    
P.show()
