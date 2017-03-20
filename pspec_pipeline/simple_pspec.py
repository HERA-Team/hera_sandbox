#! /usr/bin/env python
"""Create simple power spectrum using FFT."""

import numpy as np
import aipy as ap
import capo
import sys
import argparse
import healpy as hp
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description=('Compute Power Spectrum of files using simple FFT'
                 ' and bootstraps'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--even_files', metavar='<EVEN_FILE>', type=str, nargs='+',
                    help='List of even day files')
parser.add_argument('--odd_files', metavar='<ODD_FILE>', type=str, nargs='+',
                    help='List of odd day files')
parser.add_argument('--output', type=str, default='./',
                    help='Specifically specify out directory.')
parser.add_argument('-b', '--nboots', type=int, default=60,
                    help='Number of Bootstraps (averages) default=6')
parser.add_argument('-a', '--ant', dest='ant', default='cross',
                    help=('Select ant_pol/baselines to include. '
                          'Examples: all (all baselines) '
                          'auto (of active baselines, only i=j) '
                          'cross (only i!=j) '
                          '0,1,2 (any baseline involving listed ants) '
                          '0_2,0_3 (only listed baselines) '
                          '"(0,1)_(2,3)" (same as 0_2,0_3,1_2,2_3. '
                          'Quotes help bash deal with parentheses) '
                          '"(-0,1)_(2,-3)" (exclude 0_2,0_3,1_3 include 1_2). '
                          'Default is "cross". '
                          'Select pol by adding appropriate x or y eg 5x_6y.'))
parser.add_argument('-p', '--pol', dest='pol', default=-1,
                    help='Polarization to analyze (xx yy xy yx I)')
parser.add_argument('-c', '--chan', dest='chan', default='all', nargs='+',
                    help=('Select channels (after any delay/delay-rate '
                          'transforms) to include.  '
                          'Examples: all (all channels), '
                          '0_10 (channels from 0 to 10, including 0 and 10)'
                          '0_10_2 (channels from 0 to 10, counting by 2), '
                          '0 10 20_30 (mix of individual channels and ranges).'
                          ' Default is "all".'))
parser.add_argument('-C', '--cal', dest='cal', metavar='<CAL>.py',
                    help='Calibration file <cal>.py')
parser.add_argument('--rmbls', dest='rmbls', type=str,
                    help=('List of baselines (ex:1_4,2_33) '
                          'to remove from the power spectrum analysis.'))
parser.add_argument('--NGPS', type=int, default=5,
                    help='Number of groups for bootstrapping.')
parser.add_argument('--nside', type=int, default=512,
                    help='Nside to calculate Beam integral.')
parser.add_argument('--analytic', action='store_true',
                    help='Plot analytic noise curve.')
parser.add_argument('--Trcvr', type=float, default=200,
                    help='Receiver Temperature in Kelvin')
args = parser.parse_args()


def calc_beam(cal, nside, inttime, freq):
    """Calculate Effective Omega from Beam integral."""
    from scipy.interpolate import interp1d
    from capo import fringe
    aa1 = ap.cal.get_aa(cal, np.array([freq]))
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, pix)
    # only keep points above the horizon
    phi = phi[theta < np.pi/2.]
    theta = theta[theta < np.pi/2.]
    intpix = hp.ang2pix(nside, theta, phi)
    x, y, z = hp.pix2vec(nside, intpix)
    # need xyz in eq coords to make fringe weights
    xyz_eq = np.dot(np.linalg.inv(aa1._eq2zen), np.array([x, y, z]))
    # A is defined as the power
    A = aa1[0].bm_response((x, y, z), pol=args.pol)[0]**2

    # calculate the fringe weights for FRF'd data
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa1, (1, 4), 0, bins=bins)
    bl = aa1.get_baseline(1, 4, 'r') * freq
    fng = fringe.mk_fng(bl, xyz_eq)
    wgts = interp1d(bins, frp, kind='linear')
    fng_wgt = wgts(fng)
    fng_bm = A * fng_wgt
    return A, fng_wgt


def calc_omega(beam, nside):
    """Caclulate the Beam area and Beams^2 areas."""
    pixarea = hp.nside2pixarea(nside)  # pixel area in radianss
    Omega_p = np.sum(beam) * pixarea
    Omega_pp = np.sum(beam**2) * pixarea
    # return Omega_p**2/Omega_pp
    return Omega_p, Omega_pp


np.random.seed(0)
try:
    rmbls = []
    rmbls_list = args.rmbls.split(',')
    for bl in rmbls_list:
        i, j = bl.split('_')
        rmbls.append((int(i), int(j)))
    print 'Removing baselines:', rmbls
except:
    rmbls = []

dsets = {'even': args.even_files,
         'odd': args.odd_files}
uv = ap.miriad.UV(dsets.values()[0][0])
freqs = ap.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
nchan = uv['nchan']
sdf = uv['sdf'] * 1e9  # put in Hz
inttime = uv['inttime']
(uvw, t, ij), d = uv.read()
# load a dummy antenna array to get the conj dict
aa = ap.cal.get_aa(args.cal, np.array([.151]))
uvw = aa.get_baseline(ij[0], ij[1], src='z') * capo.cosmo_units.c * 1e-9
bl_length = np.linalg.norm(uvw)
bls, conj = capo.red.group_redundant_bls(aa.ant_layout)

try:
    frf_inttime = uv['FRF_NEBW']
except:
    frf_inttime = inttime
del(uv)


dlst = int(frf_inttime/inttime)
data_dict = {}
conj_dict = {}
flg_dict = {}
stats, lsts, data, flgs = {}, {}, {}, {}
for k in dsets.keys():
    if k not in data_dict.keys():
        data_dict[k] = {}
        flg_dict[k] = {}
    stats[k], data[k], flgs[k] = capo.miriad.read_files(dsets[k],
                                                        antstr=args.ant,
                                                        polstr=args.pol,
                                                        verbose=True)
    lsts[k] = np.array(stats[k]['lsts'])
    if rmbls:
        print "    Removing baselines:",
        for bl in rmbls:
            data[k].pop(bl, None)
            flgs[k].pop(bl, None)
            print bl,
        print '\n'
    for bl in data[k].keys():
        d = np.array(data[k][bl][args.pol])  # extract freq range
        flg = np.array(flgs[k][bl][args.pol])  # extract freq range
        conj_dict[bl] = conj[bl]
        if not conj[bl]:
            data_dict[k][bl] = np.conj(d)
        else:
            data_dict[k][bl] = d
        flg_dict[k][bl] = np.logical_not(flg)

keys = data_dict.keys()
bls_master = set()
for key in keys:  # populate list of baselines
    bls_master.update(data_dict[key].keys())
bls_master = list(bls_master)
print 'Baselines:', len(bls_master)

inds = capo.oqe.lst_align(lsts)

for k in data_dict.keys():
    for bl in data_dict[k].keys():
        data_dict[k][bl] = data_dict[k][bl][inds[k]]
        flg_dict[k][bl] = flg_dict[k][bl][inds[k]]
    lsts[k] = lsts[k][inds[k]]
lsts = lsts['even']

if np.size(args.chan) == 1:
    args.chan = args.chan[0].split(',')

for chan_range in args.chan:
    chans = ap.scripting.parse_chans(chan_range, len(freqs))
    print 'Computing Power Spectrum'
    afreqs = freqs.take(chans)
    nchan = chans.size
    fq = np.average(afreqs)
    z = capo.pspec.f2z(fq)
    aa = ap.cal.get_aa(args.cal, afreqs)
    etas = np.fft.fftshift(capo.pspec.f2eta(afreqs))  # Get delays (etas)
    kpl = etas * capo.pspec.dk_deta(z)
    beam_power, fringe_weights = calc_beam(args.cal, args.nside, inttime, fq)
    omega_p, omega_pp = calc_omega(beam_power*fringe_weights, args.nside)
    X2Y = capo.pspec.X2Y(z)/1e9  # convert 1/GHz to 1/Hz
    wavelength = capo.cosmo_units.c / (fq * 1e9)
    ubl = bl_length / wavelength
    kperp = capo.pspec.dk_du(fq) * ubl
    cnts = stats['even']['cnt'][inds['even']][:, chans]
    kmin = np.argmin(np.abs(kpl))
    kpl_fold = np.array(kpl[kmin:])
    ks = np.sqrt(kperp**2 + kpl_fold**2)
    Nt_eff = int(np.ceil(len(lsts)/dlst))

    print '\tFreq:', fq
    print '\tz:', z
    print '\tOmega_pp:', omega_pp
    print '\tOmega_p:', omega_p
    print '\tX2Y:', X2Y
    sys.stdout.flush()

    power_specs = {}
    for boot in xrange(args.nboots):
        key_set = np.copy(bls_master)
        np.random.shuffle(key_set)
        groups = np.array_split(key_set, args.NGPS, axis=0)
        groups = [[tuple(gp[np.random.choice(len(gp), replace=True)])
                  for bl in gp]
                  for gp in groups]

        data_grouped = {}
        power_grouped = {}
        for k in data_dict.keys():
            data_grouped[k] = [[data_dict[k][bl][:, chans] for bl in gp]
                               for gp in groups]
            data_grouped[k] = [np.mean(data, axis=0)
                               for data in data_grouped[k]]
            # average over groups for easie
            # FFT fo k-space
            power_grouped[k] = np.fft.ifft(np.conj(data_grouped[k]), axis=-1)
            power_grouped[k] = np.fft.ifftshift(power_grouped[k], axes=-1)

        count = 0
        power = 0
        for cnt, k1 in enumerate(data_dict.keys()):
            for k2 in data_dict.keys()[cnt+1:]:
                if k1 == k2:
                    continue
                for i in xrange(args.NGPS-1):
                    for j in xrange(i+1, args.NGPS):
                        p1 = (power_grouped[k1][i].conj()
                              * power_grouped[k2][j]).conj().real
                        # p1 = p1[::dlst]  # down select to independent modes
                        power += p1
                        count += 1
        power_specs[boot] = (power_specs.get(boot, 0)
                             + power / count * X2Y / omega_pp
                             * sdf * afreqs.size)

    averaged_pspecs, avg_fold = [], []
    for boot in xrange(args.nboots):
        # bl = np.random.choice(args.nboots, replace=True)
        # times = np.random.choice(Nt_eff, Nt_eff, replace=False)
        pos_neg = np.random.choice(1)
        # median over times
        tmp_pspec = np.median(power_specs[boot][::dlst], axis=0)
        averaged_pspecs.append(tmp_pspec)
        # split into pos and neg arrays, reverse the neg array so magnitude of
        # kpl increases
        tmp_pspec = [tmp_pspec[:kmin+1][::-1], tmp_pspec[kmin:]]
        avg_fold.append(np.mean(tmp_pspec, axis=0))  # mean over k// +/-

    pk = np.mean(averaged_pspecs, axis=0)
    pk_err = np.std(averaged_pspecs, axis=0)
    pk_fold, pk_fold_err = np.mean(avg_fold, axis=0), np.std(avg_fold, axis=0)
    k3pk = ks**3/(2*np.pi**2)*pk_fold
    k3err = ks**3/(2*np.pi**2)*pk_fold_err

    # Plot Pk and Delta^s side by side
    # 2-sigma Error bars
    plt.figure(figsize=(8, 6))
    plt.subplot(121)
    plt.errorbar(kpl, pk, 2*pk_err, fmt='k.')
    plt.grid()
    plt.ylim([-.5e7, None])

    plt.subplot(122)
    plt.errorbar(ks, k3pk, 2*k3err, fmt='k.')
    plt.grid()
    plt.yscale('log')
    plt.ylim([1e1, None])

    Tsys = 180 * (fq/.18)**-2.55 + args.Trcvr
    Tsys *= 1e3
    Npol = 2
    Nreal = 2
    folding = 2
    Nlst = len(lsts) * inttime / frf_inttime
    Nbls = len(bls_master) / np.sqrt(2) * np.sqrt(1. - 1./args.NGPS)

    # calculate the effective counts in the data, this is like ndays
    cnt_eff = 1./np.sqrt(np.ma.masked_invalid(1./cnts**2).mean())
    pk_noise = X2Y * omega_p**2/omega_pp * Tsys**2
    pk_noise /= frf_inttime * Npol * Nbls * cnt_eff
    pk_noise /= np.sqrt(Nlst * folding * Nreal)
    k3noise = ks**3/(2*np.pi**2) * pk_noise
    print '\tNbls eff:', Nbls
    print '\tNlst bins:', Nlst
    print '\tNdays eff:', cnt_eff
    print '\tPk noise [mK^2]: {0:.3e}'.format(pk_noise)
    if args.analytic:
        plt.subplot(121)
        plt.axhline(2*pk_noise, linestyle='-', color='g')
        plt.subplot(122)
        plt.plot(ks, 2*k3noise, 'g-')
    plt.suptitle('z={0:0.2f}'.format(z))
    plt.savefig('simple_pspec_z_{0:.2f}.png'.format(z), format='png')
    np.savez(args.output + 'simple_pspec_z_{0:.2f}.npz'.format(z), k=ks,
             kpl=kpl, pk=pk, pk_err=2*pk_err, k3pk=k3pk, k3err=2*k3err,
             pk_noise=pk_noise, k3noise=2*k3noise, freq=fq, redshift=z)
plt.show()
