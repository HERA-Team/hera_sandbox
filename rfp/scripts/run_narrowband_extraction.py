# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 the HERA Collaboration
# Licensed under the 2-clause BSD license.
import glob
import os
import re
import sys

from hera_qm import xrfi
from pyuvdata import UVData
from pyuvdata.utils import polstr2num

# figure out which directory the data lives in
filename = sys.argv[1]
dirname, file_basename = os.path.split(filename)
jd_pattern = re.compile("[0-9]{7}")
jd = jd_pattern.findall(file_basename)[0]
file_glob = sorted(
    glob.glob(os.path.join(dirname, "zen.{jd}.*.uvh5".format(jd=jd)))
)
file_glob = list([fname for fname in file_glob if "diff" not in fname])
if len(file_glob) == 0:
    raise FileNotFoundError(
        "Something went wrong--no files were found."
    )

# load in the data, downselect to autos and linear pols
use_pols = [polstr2num(pol) for pol in ('xx', 'yy')]
uvd = UVData()
uvd.read(file_glob, ant_str='auto', polarizations=use_pols)

# just do everything in this script; first isolate the rfi
rfi_data = np.zeros_like(uvd.data_array, dtype=uvd.data_array.dtype)

for antpairpol in uvd.get_antpairpols():
    # get indices for properly slicing through data array
    blt_inds, conj_blt_inds, pol_inds = uvd._key2inds(antpairpol)

    # approximately remove all non-rfi signal
    this_data = uvd.get_data(antpairpol)
    filt_data = xrfi.medminfilt(this_data)
    this_rfi = this_data - filt_data

    # detrend the original data to find where the stations are
    detrended_data = xrfi.detrend_medfilt(data)
    station_flags = np.where(
        detrended_data > 100, True, False
    )

    # update flags with watershed algorithm
    station_flags = xrfi._ws_flag_waterfall(
        detrended_data, station_flags, nsig=20
    )

    # update the rfi data to only keep rfi from narrowband transmitters
    # using 1 as the zero value so that it's log-friendly
    this_rfi = np.where(station_flags, this_rfi, 1)

    # update the rfi_data array appropriately
    rfi_data[blt_inds, 0, :, pol_inds[0]] = this_rfi

# TODO: update the clobber part to be reasonable
clobber = True
save_filename = os.path.join(
    dirname, "zen.{jd}.rfi_stations.uvh5".format(jd=jd)
)
uvd.data_array = rfi_data
uvd.write_uvh5(save_filename, clobber=clobber)

