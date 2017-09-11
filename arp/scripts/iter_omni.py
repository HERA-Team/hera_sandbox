#! /usr/bin/env python
from pyuvdata import UVCal, UVData
import numpy as np, pylab as plt
import hera_cal.firstcal
import aipy

fc = UVCal()
fc.read_calfits('zen.2457555.42443.xx.HH.uvc.first.calfits')
fc.convert_to_gain(run_check=False)

oc = UVCal()
oc.read_calfits('zen.2457555.42443.HH.uvc.omni.calfits')

uv = UVData()
uv.read_miriad('zen.2457555.42443.xx.HH.uvc')

d,f = hera_cal.firstcal.UVData_to_dict([uv])

aa = aipy.cal.get_aa('hsa7458_v001', np.array([.15]))
info = hera_cal.omni.aa_to_info(aa, ex_ants=[22,81])

g0_0 = {'x':{}}
for i,ai in enumerate(fc.ant_array):
    g0_0['x'][ai] = fc.gain_array[i,0,:,:,0].T
g0_0['x'][43] *= -1 # 43 appears to have a 180 deg swapped dipole

m1_0,g1_0,v1_0 = hera_cal.omni.run_omnical(d, info, g0_0)

plt.figure(0)
for i in g0_0['x']:
    plt.plot(np.angle(g1_0['x'][i][0]*np.conj(g0_0['x'][i][0])), label=i)
plt.legend()

d_cal = {}
for i,j in d:
    try: gi = g0_0['x'][i]
    except: gi = np.array([1.], dtype=np.complex64)
    try: gj = g0_0['x'][j]
    except: gj = np.array([1.], dtype=np.complex64)
    d_cal[i,j] = {'xx': (d[i,j]['xx'] / (gi * gj.conj())).astype(np.complex64)}

g_ones = {'x':{}}
for i in g0_0['x']:
    g_ones['x'][i] = np.ones_like(g0_0['x'][i])

m1_1,g1_1,v1_1 = hera_cal.omni.run_omnical(d_cal, info, g_ones)

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
