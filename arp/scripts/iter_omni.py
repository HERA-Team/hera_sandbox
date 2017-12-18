#! /usr/bin/env python
from pyuvdata import UVCal, UVData
import numpy as np, pylab as plt
import hera_cal
import aipy, sys, os

CAL = 'hsa7458_v001'
#CAL = None
ex_ants = [22,81]
switched = [43]
#POLS = ['xx','yy']
POLS = ['xx']
ubls = None

for filename in sys.argv[1:]:
    print("Reading {0}".format(filename))
    uv_in = UVData(); uv_in.read_miriad(filename)
    if uv_in.phase_type != 'drift':
        print("Setting phase type to drift")
        uv_in.unphase_to_drift()
    fqs = uv_in.freq_array[0, :] / 1e9 # Hz -> GHz
    if CAL is None: aa = hera_cal.utils.get_aa_from_uv(uv_in) # XXX ends up with ants not in HERA
    else: aa = aipy.cal.get_aa(CAL, fqs)
    print('Excluding Antennas:', ex_ants)
    if ubls != None:
        print('Using Unique Baselines:', ubls)
    # XXX determine pols from file being read
    info_fc = hera_cal.omni.aa_to_info(aa, pols=[p[0] for p in POLS], fcal=True, ubls=ubls, ex_ants=ex_ants)
    # XXX how is info_fc different from fcal=False?
    bls = [bl for bls in info_fc.get_reds() for bl in bls]
    print('Number of redundant baselines:', len(bls))

    datapack, flagpack = hera_cal.firstcal.UVData_to_dict([uv_in])
    def _apply_pi_shift(data_dict, invert_these):
        for ai, aj in data_dict.keys():
            for pol in data_dict[ai, aj].keys():
                if ((ai, pol[0]) in invert_these) ^ ((aj, pol[1]) in invert_these):
                    data_dict[ai,aj][pol] *= -1
        return data_dict
    datapack = _apply_pi_shift(datapack, switched) # XXX in place change to data
    wgtpack = {k: {p: np.logical_not(flagpack[k][p]) for p in flagpack[k]} for k in flagpack}

    # gets phase solutions per frequency.
    fc = hera_cal.firstcal.FirstCal(datapack, wgtpack, fqs, info_fc)
    sols = fc.run(finetune=True, verbose=False, average=False, window='none')
    # XXX sols has antpol object keys.  needs (ants,pol)
    def tau2gain(fqs,tau):
        fqs = np.reshape(fqs,(1,-1))
        return np.exp(2j*np.pi*fqs*tau)
    # XXX takes fqs in GHz but returns delays in s, not ns
    # XXX iteration to get offset/switched should hapeen inside fc, not outside
    g0 = {}
    g0['x'] = {antpol.ant():tau2gain(fqs,np.array(tau.T)/1e-9) for antpol,tau in sols.items() if antpol.pol() == 'x'}
    g0['y'] = {antpol.ant():tau2gain(fqs,np.array(tau.T)/1e-9) for antpol,tau in sols.items() if antpol.pol() == 'y'}
    for pol in g0.keys():
        if len(g0[pol]) == 0: del(g0[pol])
    #import IPython; IPython.embed()

    info = hera_cal.omni.aa_to_info(aa, ex_ants=ex_ants)

    #outfile = filename + 'O'
    #try: assert(not os.path.exists(outfile))
    #except(AssertionError): continue
    #uvi = UVData()
    #print 'Reading', filename
    #uvi.read_miriad(filename)
    #d,f = hera_cal.firstcal.UVData_to_dict([uvi])


    #for ai in switched: g0['x'][i] *= -1 # 43 appears to have a 180 deg swapped dipole

    m1,g1,v1 = hera_cal.omni.run_omnical(datapack, info, g0)

    if True:
        plt.figure(0)
        for i in g0['x']:
            plt.plot(np.angle(g1['x'][i][0]*np.conj(g0['x'][i][0])), label=i)
        plt.legend()
        plt.show()
    import IPython; IPython.embed()
'''

    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)

    times = []
    def mfunc(uv, p, d, f):
        _, t, (i,j) = p
        if len(times) == 0 or times[-1] != t: times.append(t)
        ti = len(times) - 1
        pi,pj = aipy.miriad.pol2str[uv['pol']]
        try: gi = g0_0[pi][i][ti] * g1_1[pi][i][ti]
        except(KeyError): gi = np.array([1.], dtype=np.complex64)
        try: gj = g0_0[pj][j][ti] * g1_1[pj][j][ti]
        except(KeyError): gj = np.array([1.], dtype=np.complex64)
        d_new = d / (gi * gj.conj())
        return p, d_new, f
        
    print 'Writing', outfile
    uvo.pipe(uvi, mfunc, raw=True)
    del(uvo)
#g0_1 = {'x':{}}
#for i in g0_0['x']: g0_1['x'][i] = g0_0['x'][i].copy()

#m1_1,g1_1,v1_1 = hera_cal.omni.run_omnical(d, info, g0_1)

plt.figure(1)
for i in g0_0['x']:
    plt.plot(np.angle(g1_1['x'][i][0]*np.conj(g_ones['x'][i][0])), label=i)
plt.legend()

plt.show()

import IPython; IPython.embed()
'''

'''
g0_2 = {'x':{}}
for i in g1_1['x']:
    phs = g1_1['x'][i][0] * np.conj(g0_1['x'][i][0])
    p = np.polyfit(fq, np.angle(phs), deg=1)
    print i, p
    plt.plot(np.angle(phs))
    plt.plot(np.polyval(p, fq))
    plt.show()
    g0_2['x'][i] = g0_1['x'][i] * np.exp(1j*np.polyval(p,fq))
m1_2,g1_2,v1_2 = hera_cal.omni.run_omnical(d,info,g0_2)
'''
