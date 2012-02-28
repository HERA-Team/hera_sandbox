#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P

filename = 'zen.2455600.70036.uvcmbRx'
uv = a.miriad.UV(filename)
fqs = n.arange(uv['nchan'])*uv['sdf'] + uv['sfreq']
del(uv)
bl_sets = ['3_11,4_12,5_13,6_14,3_19,4_20,5_21,6_22',
           '3_12,4_13,5_14,4_19,5_20,6_21',
           '3_13,4_14,5_19,6_20',
           '3_14,6_19',
           '4_11,5_12,6_13,3_20,4_21,5_22',
           '5_11,6_12,3_21,4_22',
           '6_11,3_22',
           '11_19,12_20,13_21,14_22',
           '12_19,13_20,14_21',
           '13_19,14_20',
           '11_20,12_21,13_22',
           '11_21,12_22',
           #'3_4,4_5,5_6,11_12,12_13,13_14,19_20,20_21,21_22',
           #'3_5,4_6,11_13,12_14,19_21,20_22',
           #'3_6,11_14,19_22',
]
t,d,f = C.arp.get_dict_of_uv_data(filename, ','.join(bl_sets), 'yy')
bl_sets = [[a.miriad.ij2bl(int(word.split('_')[0]),int(word.split('_')[1])) 
        for word in bllist.split(',')]
            for bllist in bl_sets]
for bl in d:
    i,j = a.miriad.bl2ij(bl)
    if j in range(16,24): d[bl] = n.conjugate(d[bl])

v = {}
for k in f: v[k] = n.logical_not(f[k])
for k in d: d[k] *= v[k]

ants, antlist = {}, []
# measured = m * param
measured, m = [], []
for bls in bl_sets:
  for cnt,bl0 in enumerate(bls):
    i0,j0 = a.miriad.bl2ij(bl0)
    # Permute indices to reflect conjugation above
    if j0 in range(16,24): i0,j0 = j0,i0
    # Build matrix indices for antennas
    for bl1 in bls[cnt+1:]:
        i1,j1 = a.miriad.bl2ij(bl1)
        # Permute indices to reflect conjugation above
        if j1 in range(16,24): i1,j1 = j1,i1
        for ant in [i0,j0,i1,j1]:
            if not ants.has_key(ant):
                antlist.append(ant)
                ants[ant] = len(ants)
                m = [line+[0.] for line in m] # add dimension to matrix along ant axis
        # Compute measured values
        tau,off,dtau,doff = 0,0,0,0
        dij01 = d[bl0] * n.conj(d[bl1])
        dij01_sum = n.sum(dij01,axis=0)
        dij01_wgt = n.sum(v[bl0]*v[bl1],axis=0)
        dij01_avg = dij01_sum / dij01_wgt
        dij01_avg_abs = n.abs(dij01_avg)
        for j in range(10):
            dij01_avg *= n.exp(-1j*(fqs*dtau+doff))
            tau += dtau; off += doff
            phs = n.where(dij01_avg_abs > 0, dij01_avg/dij01_avg_abs, 0)
            _phs = n.abs(n.fft.fft(phs))
            mx = n.argmax(_phs)
            if mx == 0:
                # Fine-tune calibration with linear fit
                valid = n.where(dij01_wgt > dij01_wgt.max()/2, 1, 0)
                fqs_val = fqs.compress(valid)
                dly = n.real(n.log(phs.compress(valid))/1j)
                dtau,doff = n.polyfit(fqs_val, dly, deg=1)
            else:
                # Pull out an integral number of phase wraps
                if mx > _phs.size/2: mx -= _phs.size
                dtau,doff = 2*n.pi*mx / (fqs[-1] - fqs[0]), 0
        off %= 2*n.pi
        print ((i0,j0),(i1,j1)),':',(tau,off)
        # Add entries to matrix equation
        measured.append(tau)
        line = [0.] * len(ants)
        # For dividing baselines, p = (pj0 - pi0) - (pj1 - pi1)
        line[ants[i0]] += -1. # pi0 = -1
        line[ants[j0]] +=  1. # pj0 =  1
        line[ants[i1]] +=  1. # pi1 =  1
        line[ants[j1]] += -1. # pj1 = -1
        m.append(line)
        # Plot stuff
        dij01_sum = n.sum(dij01,axis=0)
        dij01_wgt = n.sum(v[bl0]*v[bl1],axis=0)
        valid = n.where(dij01_wgt > dij01_wgt.max()/2, 1, 0)
        fqs_val = fqs.compress(valid)
        dly = n.real(n.log(phs.compress(valid))/1j)
        dij01_avg = dij01_sum / dij01_wgt * n.exp(-1j*(fqs*tau+off))
        dij01_avg_abs = n.abs(dij01_avg)
        phs = n.where(dij01_avg_abs > 0, dij01_avg/dij01_avg_abs, 0)
        dly = n.real(n.log(phs.compress(valid))/1j)
        P.plot(fqs_val, dly, ',', label=str(a.miriad.bl2ij(bl0)))
        P.legend(loc='best')

# Add 1 more piece of info: define first antenna to have param = 0 (arbitrary)
for i in []:
    measured.append(0.)
    m.append([0.] * len(ants))
    m[-1][i] = 1.

m = n.array(m)
measured = n.array(measured)

print ants
n.set_printoptions(threshold=n.nan)
print m
print
print measured
print
x,res,rank,s = n.linalg.lstsq(m, measured)
print x
print s

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

P.show()
