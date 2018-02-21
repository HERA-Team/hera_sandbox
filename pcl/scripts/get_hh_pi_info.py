#!/usr/bin/env python2.7
# -*- mode: python; coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import numpy as np
import itertools
import aipy
from pyuvdata import utils as uvutils
from hera_mc import cm_hookup, geo_handling
from astropy.time import Time
import astropy.constants as const

ant_nums = []
longitudes = []
latitudes = []
elevations = []

jd = 2457458
t = Time(jd, format='jd')
hookup = cm_hookup.Hookup(at_date=t)
hookup_dict = hookup.get_hookup(['HH', 'PI'])
geo = geo_handling.Handling()

# get cofa info
cofa_loc = geo.cofa()[0]
latitude = cofa_loc.lat * np.pi / 180.
longitude = cofa_loc.lon * np.pi / 180.
altitude = cofa_loc.elevation
location = latitude, longitude, altitude

# get antennas and their locations (lat, lon, alt)
for hkey in hookup_dict['hookup'].keys():
    for pkey in hookup_dict['hookup'][hkey]:
        ant_number = int(hookup_dict['hookup'][hkey][pkey][0].downstream_part[1:])
        if ant_number not in ant_nums:
            stn_name = hookup_dict['parts_epoch']['path'][0]
            stn = cm_hookup.get_parts_from_hookup(stn_name, hookup_dict)[hkey][pkey]
            fnd_list = geo.get_location([stn[0]], t, station_types=geo.station_types)
            fnd = fnd_list[0]
            lat = fnd.lat
            lon = fnd.lon
            alt = fnd.elevation

            ant_nums.append(ant_number)
            longitudes.append(lon)
            latitudes.append(lat)
            elevations.append(alt)

# convert to ecef (in meters)
ecef_positions = uvutils.XYZ_from_LatLonAlt(np.array(latitudes) * np.pi / 180.,
                                            np.array(longitudes) * np.pi / 180.,
                                            elevations)
rotecef_positions = uvutils.rotECEF_from_ECEF(ecef_positions.T,
                                              cofa_loc.lon * np.pi / 180.)

# make an array of antenna positions in the form miriad is expecting
c_ns = const.c.to('m/ns').value
nants = max(ant_nums) + 1
antpos = np.zeros([nants, 3])
for i, ant in enumerate(ant_nums):
    antpos[ant, :] = rotecef_positions[i, :] / c_ns

# miriad magic numbers
n_chans = 1024
sdf = 9.765625e-05
sfreq = 0.1 # bandwidth in GHz

# make an aa object
freqs = np.arange(n_chans, dtype=np.float) * sdf + sfreq
beam = aipy.phs.Beam(freqs)
ants = [aipy.phs.Antenna(a[0], a[1], a[2], beam) for a in antpos]
aa = aipy.phs.AntennaArray[ants=ants, location=location)

# loop over miriad files

# XXX: DEFINE LIST OF MIRIAD FILES
for fn in list_of_miriad_files:
    uvd = UVData()
    uvd.read_miriad(fn)

    # set the telescope location
    uvd.telescope_location_lat_lon_alt = np.array(cofa_loc.lat * np.pi / 180.,
                                                  cofa_loc.lon * np.pi / 180.,
                                                  cofa_loc.elevation)

    # loop over aa object
    idx = 0
    antpos = np.zeros((len(aa), 3))
    ants_telescope = []
    for iant, ant in enumerate(aa):
        # test to see if antenna is "far" from center of the Earth
        if np.linalg.norm(ant.pos) > 1e6:
            # convert from ns -> m
            pos = ant.pos * c_ns

            # rotate from rotECEF -> ECEF
            pos_ecef = uvutils.ECEF_from_rotECEF(pos, longitude)

            # subtract off array location, to get just relative positions
            rel_pos = pos_ecef - uv.telescope_location

            # save in array
            antpos[idx, :] = rel_pos

            # also save antenna number to list of antennas
            ants_telescope.append(iant)

            # increment counter
            idx += 1

    # set the antenna information
    uvd.Nants_telescope = len(ants_telescope)
    uvd.antenna_numbers = np.asarray(ant_nums)
    uvd.antenna_positions = np.array(antpos[:idx, :])

    # generate uvw
    uvw = []
    # get first frequency in aa object, and convert from GHz -> Hz
    try:
        freq = aa.get_freqs()[0] * 1e9
    except IndexError:
        freq = aa.get_freqs() * 1e9
    # get wavelength in meters
    lamb = const.c.to('m/s').value / freq
    # get baselines
    bls = np.asarray(itertools.product(ant_nums, ant_nums))
    bls = sorted(map(uv.antnums_to_baseline, bls[:, 0], bls[:, 1]))
    for t in range(uv.Ntimes):
        for bl in bls:
            uvw.append(aa.gen_uvw(
                *uv.baseline_to_antnums(bl), src='z').reshape(3, -1))
    # multiply by wavelength
    uv.uvw_array = np.array(uvw).reshape(-1, 3) * lamb

    # DEFINE WHAT NEW FILE NAME SHOULD BE
    uvd.write_miriad('new_file.uv')
