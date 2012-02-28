#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P

filename = 'zen.2455600.70036.uvcmbRx'
uv = a.miriad.UV(filename)
fqs = n.arange(uv['nchan'])*uv['sdf'] + uv['sfreq']
del(uv)
INT = 113
CHAN = 605
WIERINGA = False
USE_VERTICAL = True
if not WIERINGA:
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
    ]
    if USE_VERTICAL:
        bl_sets += [
               '3_4,4_5,5_6,11_12,12_13,13_14,19_20,20_21,21_22',
               '3_5,4_6,11_13,12_14,19_21,20_22',
               '3_6,11_14,19_22',
        ]
else:
    bl_sets = ['1_2,2_3,3_4,4_5',
               '1_3,2_4,3_5',
               '1_4,2_5',
    ]
#t,d,f = C.arp.get_dict_of_uv_data(filename, ','.join(bl_sets), 'yy')
t,d,f = [],{},{}
bl_sets = [[a.miriad.ij2bl(int(word.split('_')[0]),int(word.split('_')[1])) 
        for word in bllist.split(',')]
            for bllist in bl_sets]
for bl in d:
    i,j = a.miriad.bl2ij(bl)
    if j in range(16,24) and not i in range(16,24): d[bl] = n.conjugate(d[bl])

v = {}
for k in f: v[k] = n.logical_not(f[k])
for k in d: d[k] *= v[k]

# Build index of parameters to solve for
# meas = m * prms
prms, prmlist = {}, []
for cnt, bls in enumerate(bl_sets):
    # Add type of visibility as a parameter
    prms['vis%d'%cnt] = len(prmlist)
    prmlist.append(1+0j)
    for bl in bls:
        i,j = a.miriad.bl2ij(bl)
        # Permute indices to reflect conjugation above
        if j in range(16,24) and not i in range(16,24): i,j = j,i
        # Add antenna parameters
        for ant in [i,j]:
            if not prms.has_key('ant%d'%ant):
                prms['ant%d'%ant] = len(prms)
                prmlist.append(1+0j)
        
# Make fake data
if not WIERINGA:
    fake = {
    'ant3':0.0,
    'ant4':0.4,
    'ant5':0.5,
    'ant6':0.6,
    'ant11':0.11,
    'ant12':0.12,
    'ant13':0.13,
    'ant14':0.14,
    'ant19':0.19,
    'ant20':0.20,
    'ant21':0.21,
    'ant22':0.22,
    'vis0': 1.0,
    'vis1': 0.9,
    'vis2': 0.8,
    'vis3': 0.7,
    'vis4': 1.1,
    'vis5': 1.2,
    'vis6': 1.3,
    'vis7': 2,
    'vis8': 1.9,
    'vis9': 1.8,
    'vis10': 2.1,
    'vis11': 2.2,
    }
    if USE_VERTICAL:
        fake.update({
            'vis12':-0.1,
            'vis13':-0.2,
            'vis14':-0.3,
        })
else:
    fake = { # Weiringa 
        'ant1':-0.2,
        'ant2':-0.1,
        'ant3': 0.0,
        'ant4': 0.1,
        'ant5': 0.2,
        'vis0': 1.0,
        'vis1': 2.0,
        'vis2': 3.0,
    }

# Build matrix relating measurements to parameters
meas, m = [], []
for cnt, bls in enumerate(bl_sets):
    for bl in bls:
        i,j = a.miriad.bl2ij(bl)
        # Permute indices to reflect conjugation above
        if j in range(16,24) and not i in range(16,24): i,j = j,i
        #print n.abs(d[bl][INT,CHAN])
        #vis_ij = n.log(d[bl][INT,CHAN]/n.abs(d[bl][INT,CHAN])) / 1j
        vis_ij = fake['ant%d'%j] - fake['ant%d'%i] + fake['vis%d'%cnt] # Difference parameters
        #vis_ij = fake['ant%d'%j] + fake['ant%d'%i] + fake['vis%d'%cnt] # Additive parameters
        print (i,j), '= %f - %f + %f' % (fake['ant%d'%j], fake['ant%d'%i], fake['vis%d'%cnt]), '=', vis_ij
        mline = n.zeros(len(prmlist), dtype=n.float)
        mline[prms['ant%d'%j]] = 1.
        mline[prms['ant%d'%i]] = -1. # for prms that are conjugated
        #mline[prms['ant%d'%i]] = 1. # for prms that are not conjugated
        mline[prms['vis%d'%cnt]] = 1.
        meas.append(vis_ij.real)
        m.append(mline)
# Additional info: absolute (arbitrary) phase reference
mline = n.zeros(len(prmlist), dtype=n.float)
mline[prms['ant3']] = 1
meas.append(0.)
m.append(mline)

# Additional information: absolute phase reference for sky (NS and EW)
mline = n.zeros(len(prmlist), dtype=n.float)
mline[prms['vis0']] = 1
meas.append(1.)
m.append(mline)
mline = n.zeros(len(prmlist), dtype=n.float)
mline[prms['vis12']] = 1
meas.append(-.1)
m.append(mline)

#print prms
m = n.array(m)
actuals = n.zeros(len(prmlist), dtype=n.float)
for k in prms:
    actuals[prms[k]] = fake[k]
print actuals
meas = n.array(meas)
print m.shape
n.set_printoptions(threshold=n.nan)
print m
#print m
print
print n.around(meas, 2)
print n.around(n.dot(m, actuals), 2)
#print
x,res,rank,s = n.linalg.lstsq(m, meas)
print n.around(n.dot(m,x), 2)
print
print x
print s
print '-'*60
for k in prms: print k, n.around(x[prms[k]],2), actuals[prms[k]]
'''
m = [line+[0.] for line in m] # add dimension to matrix along ant axis
        
    for bl1 in bls[cnt+1:]:
        i1,j1 = a.miriad.bl2ij(bl1)
        # Permute indices to reflect conjugation above
        if j1 in range(16,24): i1,j1 = j1,i1
        for ant in [i0,j0,i1,j1]:
            if not ants.has_key(ant):
                antlist.append(ant)
                ants[ant] = len(ants)
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
'''
