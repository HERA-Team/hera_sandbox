#! /usr/bin/env python
"""Compute the redudnant visibilities from a uvfits file."""

import numpy as np
from pyuvdata import UVData, UVCal
import hera_cal
import argparse
import sys

parser = argparse.ArgumentParser(
    description='Compute the redundant Vilisbilites and save to new file.')
parser.add_argument('files', metavar='<FILE>', type=str, nargs='*')
parser.add_argument('--noise_files', metavar='<NOISE>', type=str, nargs='*')
args = parser.parse_args()


def rms(x, **kwargs):
    """Compute RMS for complex values along given axis."""
    return np.sqrt(np.mean(x.conj() * x, **kwargs))


def repeat_visibilities(uvdata_obj, noise_array=False):
    """Replicate the visibilites for each redundant group.

    Setting noise_array=True will compute the RMS for each time step
    and compute new noise realization.
    """
    # Get a list of redundant base groups and the ant_pairs in each group
    aa_object = hera_cal.utils.get_aa_from_uv(uvdata_obj)
    aa_info = hera_cal.omni.aa_to_info(aa_object)
    reds = aa_info.get_reds()
    # Find the ant pairs stored in our object
    # and the redundant group the belong to.
    uv_ij = [tuple(ij) for ij in
             np.unique(zip(uvdata_obj.ant_1_array, uvdata_obj.ant_2_array),
             axis=0)]
    flat_red = [bl for group in reds for bl in group]
    # Use all the non_redundant baselines to
    # form the base object we will return
    non_red = [tuple(ij) for ij in uv_ij
               if np.logical_and(ij not in flat_red, ij[::-1] not in flat_red)]
    full_observation = uvdata_obj.select(ant_pairs_nums=non_red, inplace=False)
    # These parameters should be the same for all redundant groups
    nfreqs = full_observation.Nfreqs
    ntimes = full_observation.Ntimes

    for cnt, ij in enumerate(uv_ij):
        for j in xrange(len(reds)):
            if not np.logical_or(ij in reds[j], ij[::-1] in reds[j]):
                continue
            else:
                red_group = uvdata_obj.select(ant_pairs_nums=ij, inplace=False)
                # Number of baselines in the Redundant Group
                nbls = len(reds[j])

                if noise_array:
                    # Create a new noise realization for
                    # each frequency as a function of time
                    noise_rms = rms(red_group.data_array,
                                    axis=-2, keepdims=True)
                    # Each baseline should have the same RMS for each time step
                    noise_rms = np.repeat(noise_rms, nbls, axis=0)
                    shape = np.asarray(noise_rms.shape)
                    shape[-2] = nfreqs
                    # Sqrt(2) used to preserve total
                    # amplitude of noise across real and imaginary
                    noise_realization = (np.random.normal(size=shape) +
                                         1j * np.random.normal(size=shape))
                    red_group.data_array = noise_rms / np.sqrt(2) \
                        * np.array(noise_realization)
                else:
                    red_group.data_array = np.repeat(red_group.data_array,
                                                     nbls, axis=0)
                # Replicate all the necessary data and meta-data for
                # the observation
                red_group.flag_array = np.repeat(red_group.flag_array,
                                                 nbls, axis=0)
                red_group.nsample_array = np.repeat(red_group.nsample_array,
                                                    nbls, axis=0)
                red_group.uvw_array = np.repeat(red_group.uvw_array,
                                                nbls, axis=0)
                red_group.lst_array = np.repeat(red_group.lst_array,
                                                nbls, axis=0)
                red_group.time_array = np.repeat(red_group.time_array,
                                                 nbls, axis=0)
                ant_num_array = np.array([np.asarray(reds[j])[:, 0],
                                          np.asarray(reds[j])[:, 1]])

                ant_num_array = np.tile(ant_num_array, (1, ntimes))
                baseline_array = red_group.antnums_to_baseline(*ant_num_array)
                red_group.baseline_array = baseline_array
                red_group.ant_1_array = ant_num_array[0]
                red_group.ant_2_array = ant_num_array[1]
                red_group.Nants_data = len(np.unique(ant_num_array))
                red_group.Nblts = nbls * ntimes
                red_group.Nbls = nbls
                full_observation.__add__(red_group, inplace=True,
                                         check_extra=True,
                                         run_check_acceptability=True,
                                         run_check=True)
    return full_observation

if not args.files and not args.noise_files:
    print "No Input Files given."""
    sys.exit()

if args.files:
    print 'Processing File: ',
    for filename in args.files:
        print filename,
        data = UVData()
        try:
            data.read_uvfits(filename)
        except IOError:
            filename += '.uvfits'
            data.read_uvfits(filename)

        full_data = repeat_visibilities(data, noise_array=False)
        outname = filename.split('.')
        outname[-2] += '_expanded'
        outname = '.'.join(outname)
        full_data.write_uvfits(outname)
    print 'Done'

if args.noise_files:
    print 'Processing Noise File: ',
    for filename in args.noise_files:
        noise = UVData()
        try:
            noise.read_uvfits(filename)
        except IOError:
            filename += '.uvfits'
            noise.read_uvfits(filename)
        print filename,
        full_noise = repeat_visibilities(noise, noise_array=True)
        outname = filename.split('.')

        outname[-2] += '_expanded'
        outname = '.'.join(outname)
        full_noise.write_uvfits(outname)
    print 'Done'
