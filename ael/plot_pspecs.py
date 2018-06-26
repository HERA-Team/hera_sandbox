# coding: utf-8
import hera_pspec as hp
import numpy as np
import matplotlib.pyplot as plt
from eorsky import comoving_voxel_volume

saved = hp.PSpecContainer('hera_3bl_fhd_red.hdf5')
#ps = saved.get_pspec("eor_flatsig_mwabeam")[0]   # 0th spw
ps = saved.get_pspec("eor_bubble_nside512_sig2")[0]   # 0th spw
blpairs = [ps.blpair_to_antnums(p) for p in np.unique(ps.blpair_array)]

avg_ps = ps.average_spectra(blpair_groups=[blpairs], inplace=False, time_avg=True)
kpara = avg_ps.get_kparas(0)
key = (0, blpairs[0], 'pI')
power = np.abs(avg_ps.get_data(key)).T
#import IPython; IPython.embed()

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.plot(kpara, power)
#ax.set_yscale('log')
#ax = hp.plot.delay_spectrum(ps, [blpairs,], spw=0, pol='pI', average_blpairs=True, average_times=True, delay=False)

ref_level = 2193384.8800389073 * 0.5

nside=512.
dnu = (ps.beam_freqs[1] - ps.beam_freqs[0]) / 1e6
om = 4*np.pi/(12.*nside**2)
Z = np.mean(1420e6/ps.beam_freqs - 1)
ref_level = comoving_voxel_volume(Z, dnu, om) * 4.0e6

onehour = 1/24.
ps_avgs = []
ps_vars = []
times = []
frac=3.
for hi in np.arange(1,frac*24):
    inds = np.where(ps.time_avg_array < (ps.time_avg_array[0] + hi*onehour/frac))
    times.append(hi/frac)
    print len(inds[0])
    ps_avgs.append(np.mean(ps.data_array[0][inds[0],:,0], axis=0))
    ps_vars.append(np.var(ps.data_array[0][inds[0],:,0], axis=0))

#ax.axhline(y=ref_level)
#plt.xlabel(r"k$_\parallel$ Mpc$^{-1}$")
#plt.ylabel(r"P(k) mK$^2$ Mpc$^3$")
#plt.show()
