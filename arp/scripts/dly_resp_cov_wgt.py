#! /usr/bin/env python
import numpy as np, pylab as plt
import capo, aipy, scipy.interpolate, random
import capo.oqe as oqe
import sys

HERA_BM_EFF_150 = 0.0905 # at 150 MHz, from beam_area.py X4Y2_490_150.hmap for \Omega_eff
HERA_BM_NOS_150 = 0.0401 # at 150 MHz, from beam_area.py X4Y2_490_150.hmap for \Omega_eff

#CH0,NCHAN = 30, 61
CH0,NCHAN = 98, 61
#CH0,NCHAN = 1, 253
#CH0,NCHAN = 64, 127

def jy2T(fq, bm=HERA_BM_NOS_150):
    '''Return [mK] / [Jy] for a beam size (in sr) and frequency (in GHz).'''
    lam = aipy.const.c / (fq * 1e9)
    return 1e-23 * lam**2 / (2 * aipy.const.k * bm) * 1e3

def gen_eor_mdl(ks, Pk, shape, fq_axis=1, fq0=.1, fq1=.2, bm_eff=HERA_BM_EFF_150, bm_nos=HERA_BM_NOS_150, kpr=.06):
    fqs = np.linspace(fq0,fq1,shape[fq_axis]) # GHz
    B = fq1 - fq0 # GHz
    z = capo.pspec.f2z(np.average(fqs))
    scalar = capo.pspec.X2Y(z) * bm_eff * B
    mdl = scipy.interpolate.interp1d(ks, Pk, bounds_error=False, fill_value=0)
    etas = np.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
    kpls = etas * capo.pspec.dk_deta(z)
    ks = np.sqrt(kpls**2 + kpr**2)
    Pk = mdl(ks)
    dspec_amp = np.sqrt(Pk / scalar)
    sh = [1] * len(shape); sh[fq_axis] = -1
    dspec_amp.shape = sh
    jy2mK_HERA = jy2T(fqs, bm=bm_nos); jy2mK_HERA.shape = sh
    vis = oqe.noise(shape).astype(np.complex128) * dspec_amp
    vis_mK = np.fft.fft(vis, axis=fq_axis)
    vis_Jy = vis_mK / jy2mK_HERA
    return vis_Jy

def gen_constraint(tau, achr, eor, tcut=200.):
    mn = {}
    for i,ti in enumerate(tau):
        for j,tj in enumerate(tau):
            if abs(tj) < tcut: continue
            dt = np.abs(np.around(5*(tj-ti), 0)/5)
            #mn[dt] = min(abs(eor[j] / achr[i]), mn.get(dt,1)) # if in amplitude
            mn[dt] = min(eor[j] - achr[i], mn.get(dt,0)) # if in dB
    dts = mn.keys(); dts.sort()
    mns = np.array([mn[dt] for dt in dts])
    return np.array(dts), mns

#def dB(sig): return 10*np.log10(np.average(sig.real, axis=1)) / 2. # from mK^2 to mK
def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1))) / 2. # from mK^2 to mK
#def dB(sig): return 10*np.log10(np.abs(np.average(np.abs(sig), axis=1))) / 2. # from mK^2 to mK
#def dB(sig): return 10*np.log10(np.abs(np.median(sig, axis=1))) / 2. # from mK^2 to mK

def wideband_clean(fg):
    return fg
    window = aipy.dsp.gen_window(fg.shape[-1], 'blackman-harris')
    window.shape = (1,-1)
    _fg = np.fft.ifft(window*fg, axis=1)
    #_fg[:,:8] = 0; _fg[:,-7:] = 0
    #_fg[:,:15] = 0; _fg[:,-14:] = 0
    _fg[:,:21] = 0; _fg[:,-20:] = 0
    #_fg[:,:5] = 0; _fg[:,-4:] = 0
    return np.fft.fft(_fg,axis=1)/window

'''
'bl': 3x3 baseline vectors where each row is a baseline vector in meters
'tau': 512-element array of delays in seconds
'lst': 80-element array of LST in hours
'skyvis': 3x512x80 Blackman-Harris window response squared applied on delay spectra with no foreground cleaning or removal. Units of P(k) K^2 (Mpc/h)^3
'cc_skyvis_net': 3x512x80 Blackman-Harris window response squared applied on delay spectra and cleaned with residuals added back. Units of P(k) K^2 (Mpc/h)^3
'''

npz = np.load(sys.argv[-1])
dk_deta = capo.pspec.dk_deta(9.)
bls = npz['bl']
bl0 = 0 # do shortest baselines
#bl0 = 3 # do 2nd shortest baselines
print bls[bl0:bl0+3]
WINDOW = 'blackman-harris'
try:
    fqs = npz['freq'] / 1e9 # GHz
    sky = npz['skyvis_freq'].astype(np.complex128) * 1e3 * 3e3 # mK (bl,fq,t)
    sig = npz['vis_noise_freq'].astype(np.complex128) * 1e3 * 1e2 # mK (bl,fq,t)
    window = aipy.dsp.gen_window(128, WINDOW); window.shape = (1,-1,1)
except(KeyError):
    fqs = npz['f'] / 1e9 # GHz
    sky = npz['skyvis'].astype(np.complex128) # Jy (bl,fq,t)
    sky = wideband_clean(sky)
    npz = np.load('21cmfastmodel_150MHz.npz')
    #npz = np.load('lidzmodel_150MHz.npz')
    sig = gen_eor_mdl(npz['k'], npz['Pk'].astype(np.complex128) * 1e6, sky.shape)
    sig = wideband_clean(sig)
    window = aipy.dsp.gen_window(256, WINDOW); window.shape = (1,-1,1)

#pW = 1.6*2*np.abs(np.fft.fftshift(np.fft.ifft(window*sky, axis=1), axes=1))**2
plt.subplot(211); capo.plot.waterfall(sky[0], drng=6); plt.colorbar()
plt.subplot(212); capo.plot.waterfall(sig[0], drng=3); plt.colorbar()
#plt.subplot(212); capo.plot.waterfall(sky2[0], drng=8); plt.colorbar()
plt.show()

#npz_achr = np.load('fgdps_achromatic.npz')
#sky2 = npz_achr['skyvis'].astype(np.complex128) * 1e6 # mK^2 (bl,dly,t)


npz = np.load('delayspectrum.npz')
tau_meas = npz['delay']
ref_meas = npz['amp']; ref_meas -= ref_meas[np.where(tau_meas == 0)]


#fqs = np.linspace(.1,.2,NCHAN)
#fqs = np.linspace(.130,.160,NCHAN)
fqs = fqs[CH0:CH0+NCHAN]
dly = np.fft.fftfreq(NCHAN, fqs[1]-fqs[0])
k0 = ('even',(0,1),'I')
k1 = ('even',(1,2),'I')
k2 = ('even',(2,3),'I')
k3 = ('even',(3,4),'I')
fg,eor = {},{}
if True:
    fg[k0]  = sky[bl0+0,CH0:CH0+NCHAN].T
    fg[k1]  = sky[bl0+1,CH0:CH0+NCHAN].T
    fg[k2]  = sky[bl0+2,CH0:CH0+NCHAN].T
    eor[k0] = sig[bl0+0,CH0:CH0+NCHAN].T
    eor[k1] = sig[bl0+1,CH0:CH0+NCHAN].T
    eor[k2] = sig[bl0+2,CH0:CH0+NCHAN].T
else:
    NSAMP = 100
    ts = np.linspace(0,2*np.pi,NSAMP)
    NFG = 400
    scalar = 1

    def mk_dly_profile(dbs, slopes):
        prf = 0.
        for db,slope in zip(dbs,slopes):
            prf += 10**((db + slope * np.abs(dly))/10.)
        phs = np.random.uniform(0,2*np.pi,size=NCHAN)
        phs[0] = 0; phs[1:NCHAN/2+1] = phs[NCHAN/2+1:]
        return np.fft.fft(prf*np.exp(1j*phs))
    #bp = mk_dly_profile([0,-10],[-1.,-.1]))
    fg = 0
    for cnt in xrange(NFG):
        fg_ch = 1e3 * (fqs/.150)**np.random.uniform(-2,0); fg_ch.shape = (-1,1)
        bm = mk_dly_profile([np.random.uniform(-10,0),np.random.uniform(-40,0)],[np.random.uniform(-1,-2),np.random.uniform(-.1,-.2)]); bm.shape = (-1,1)
        fg_t = np.sin((cnt+1)*ts); fg_t.shape = (1,-1)
        fg += bm * fg_ch * fg_t
    eor = .01*oqe.noise(size=(NCHAN,NSAMP))

#f = 0.3 # inject less eor so pspec goes below the eor signal used for computing ratios
f = 1 # inject less eor so pspec goes below the eor signal used for computing ratios
dat = {}
for k in fg: dat[k] = (fg[k] + f*eor[k])
NSAMP = fg[k].shape[0]
dat_cut = {}
for k in dat: dat_cut[k] = np.concatenate([dat[k][:54],dat[k][65:]], axis=0)
ds = oqe.DataSet(dsets=dat)
#ds = oqe.DataSet(dsets=dat_cut)

prf_c_final = {.1:100, .15:100, .2:100}
prf_w_final = {.1:100, .15:100, .2:100}
for boot in xrange(1):
    print boot
    # gather C,iC matrices for baselines to try cross-application
    ds.clear_cache()
    ds.set_data(dat)
    Cs,iCs = {},{}
    for k in dat:
        #Cs[k] = ds.C(k)
        #Cs[k] = sum([ds.C(ki)+0*np.identity(NCHAN) for ki in dat])
        Cs[k] = sum([ds.C(ki)+3e-6*np.identity(NCHAN) for ki in dat if ki != k])
        #Cs[k] = sum([ds.C(ki)+1e-6*np.identity(NCHAN) for ki in dat if ki != k])
        #Cs[k] = sum([ds.C(ki)+1e-4*np.identity(NCHAN) for ki in dat if ki != k])
        #Cs[k] = sum([ds.C(ki)+0*np.identity(NCHAN) for ki in dat if ki != k])
        #ds.set_C({k:Cs[k]+1e-1*np.identity(NCHAN)}) # regularize a bit with some diagonal
        ds.set_C({k:Cs[k]})
        iCs[k] = ds.iC(k)

    ds.set_data(dat_cut)
    #ds.set_data(dat)
    tau = np.fft.fftshift(dly)
    win1 = aipy.dsp.gen_window(NCHAN, 'blackman-harris'); win1.shape = (-1,1)
    win2 = aipy.dsp.gen_window(NCHAN, 'blackman-harris')**1.5; win2.shape = (-1,1)
    for ki in iCs:
        qI = ds.q_hat(ki,ki,use_cov=False)
        FI = ds.get_F(ki,ki,use_cov=False)
        MI,WI = ds.get_MW(FI, mode='I')
        pI = ds.p_hat(MI,qI)
        ds.set_data(eor)
        qI_eor = ds.q_hat(ki,ki,use_cov=False)
        pI_eor = ds.p_hat(MI,qI_eor)
        ds.set_data(dat_cut)
        pW1 = 1.6*2*np.abs(np.fft.fftshift(np.fft.ifft(win1*dat_cut[ki].T, axis=0), axes=0))**2
        pW2 = 2.4*2*np.abs(np.fft.fftshift(np.fft.ifft(win2*dat_cut[ki].T, axis=0), axes=0))**2
        plt.figure(1)
        plt.plot(tau, dB(pI), 'b', label='I')
        plt.plot(tau, dB(pW1), 'g', label='W')
        plt.plot(tau, dB(pW2), 'g', label='W')
        plt.plot(tau, dB(pI_eor), 'k', label='E')
        ds.set_iC({ki:iCs[ki]})
        qC = ds.q_hat(ki,ki)
        FC = ds.get_F(ki,ki)
        #MC,WC = ds.get_MW(FC, mode='I')
        MC,WC = ds.get_MW(FC, mode='F^-1/2')
        #MC,WC = ds.get_MW(FC, mode='F^-1')
        pC = ds.p_hat(MC,qC)
        plt.figure(1)
        plt.plot(tau, dB(pC), 'r', label='C')
        #print k_1/dk_deta
        #print k_2/dk_deta
        #plt.plot(k_1/dk_deta, 10*np.log10(eor1))
        #plt.plot(k_2/dk_deta, 10*np.log10(eor2))
        plt.figure(2)
        for kcut in [.1,.15,.2]:
            tcut = kcut / dk_deta
            tau_c,prf1_c = gen_constraint(tau, np.where(dB(pW1) > -27, dB(pW1), -50), dB(pI_eor), tcut=tcut)
            prf_w_final[kcut] = np.where(prf1_c < prf_w_final[kcut], prf1_c, prf_w_final[kcut])
            tau_c,prf2_c = gen_constraint(tau, dB(pW2), dB(pI_eor), tcut=tcut)
            plt.plot(tau_c, prf1_c, 'g', label='W-%0.2f' % (kcut))
            plt.plot(tau_c, prf2_c, 'g', label='W-%0.2f' % (kcut))
            tau_c, prf_c = gen_constraint(tau, np.where(dB(pC) > -35, dB(pC), -50), dB(pI_eor), tcut=tcut)
            prf_c_final[kcut] = np.where(prf_c < prf_c_final[kcut], prf_c, prf_c_final[kcut])
            plt.plot(tau_c, prf_c, 'r', label='C-%0.2f' % (kcut))

for kcut in prf_c_final:
    plt.plot(tau_c, prf_c_final[kcut], 'r', label='C-%0.2f' % (kcut))
    plt.plot(tau_c, prf_w_final[kcut], 'g', label='C-%0.2f' % (kcut))
#plt.figure(1); plt.legend(loc='best')
plt.figure(2)
plt.plot(tau_meas, ref_meas, 'k')
#plt.legend(loc='best')
plt.xlim(0,500)
plt.show()

sav = {'tau_cstr':tau_c, 'tau_meas':tau_meas, 'refl_meas':ref_meas}
for kcut in prf_c_final:
    sav['refl_cstr_dspec_k%0.2f' % kcut] = prf_w_final[kcut]
    sav['refl_cstr_oqe_k%0.2f' % kcut] = prf_c_final[kcut]

print 'Writing refl_cstr.npz'
np.savez('refl_cstr.npz', **sav)
