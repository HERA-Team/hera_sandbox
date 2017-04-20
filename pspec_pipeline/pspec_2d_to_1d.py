#! /usr/bin/env python
"""Transforms injects to power spectrum.

Take as input a directory pointing to the output from sim_sigloss
output a power spectrum for each injection
"""
import sys
from glob import glob
import argparse
from capo.eor_results import read_bootstraps_dcj, average_bootstraps
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np

parser = argparse.ArgumentParser(
            description='Calculate power spectra for a run from sigloss_sim')
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of files to average')
parser.add_argument('--t_eff', type=int, required=True,
                    help=('effective length of integration timescale'
                          ' in # of integrations (equals FRF width/inttime)'))
#parser.add_argument('--bl_length', type=float,  required=True,
#                    help='length of baseline in meters')
parser.add_argument('--sub_pCv', action='store_true',
                    help='Remove pCv back from pC before averaging')
parser.add_argument('--outfile', type=str, default='',
                    help='Specifically specify out directory.')
parser.add_argument('--nboots', type=int, default=100,
                    help='Number of Bootstraps (averages) default=100')
# parser.add_argument('--mode', dest='mode', choices=['prob','excess'],
#                     default='prob',
#                     help='limit estimation method')
args = parser.parse_args(sys.argv[1:])

pspecs = read_bootstraps_dcj(args.files)
Nlstbins = np.shape(pspecs['pC'])[-1]
# get the number of lst integrations in the dataset
Neff_lst = np.ceil(Nlstbins/args.t_eff)
# compute the effective number of LST bins
# print Neff_lst
#  lets round up because this 'N' is only approximate

for key in pspecs.keys():
    if key == 'cmd':
        continue
    try:
        pspecs[key] = np.ma.masked_invalid(np.array(pspecs[key]))
    except:
        import ipdb
        ipdb.set_trace()
if args.sub_pCv:
    pspecs['pC'] -= pspecs['pCv']

# compute Pk vs kpl vs bootstraps
pk_pspecs = average_bootstraps(pspecs, Nt_eff=Neff_lst,
                               Nboots=args.nboots, avg_func=np.mean)

# Compute |k|

wavelength = cosmo_units.c/(pspecs['freq']*1e9)
ubl = args.bl_length/wavelength
kperp = dk_du(pspecs['freq'])*ubl
print "freq = ", pspecs['freq']
print "kperp = ", kperp
pk_pspecs['k'] = np.sqrt(kperp**2 + pk_pspecs['kpl_fold']**2)

# apply corrections to all the various channels
pspec_channels = ['pC', 'pI', 'pCv', 'pIv']
corrections = [1/np.sqrt(2),  # the median overestimates by sqrt(2)
               1.39]        # beam^2
# for chan in pspec_channels:
#     for c in [chan,chan+'_fold']: #no need to correct error bars
#         for correction in corrections:
#             pk_pspecs[c] *= correction
# make a pspec file
for key in pk_pspecs.keys():
    if key == 'cmd':
        continue
    try:
        pk_pspecs[key].fill_value = 0
        pk_pspecs[key] = pk_pspecs[key].filled()
    except:
        import ipdb
        ipdb.set_trace()
print injection_dir+'/pspec_pk_k3pk.npz'
np.savez(injection_dir+'/pspec_pk_k3pk.npz', **pk_pspecs)
