# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 the HERA Collaboration
# Licensed under the 2-clause BSD license.
import os
import re
import sys

import numpy as np
from hera_qm import xrfi
from pyuvdata import UVData

# Some file setup stuff.
auto_file = sys.argv[1]
save_dir = sys.argv[2]
pattern = re.compile("[0-9]{7}")
jd = pattern.findall(auto_file)[0]
uvd = UVData()
uvd.read(auto_file)

# Do the detrending and flagging.
data_shape = uvd.data_array.shape
detrended_rfi_data = np.zeros(data_shape, dtype=np.float)
rfi_flags = np.zeros(data_shape, dtype=np.bool)
medfilt_nsig = 40
ws_nsig = 5

for antpairpol in uvd.get_antpairpols():
    # get indices for properly slicing through data array
    blt_inds, conj_blt_inds, pol_inds = uvd._key2inds(antpairpol)
    if len(blt_inds) > 0:
        this_slice = (blt_inds, 0, slice(None), pol_inds[0])
    else:
        this_slice = (conj_blt_inds, 0, slice(None), pol_inds[1])

    # detrend the original data to find where the stations are
    this_data = uvd.get_data(antpairpol).real
    detrended_data = xrfi.detrend_medfilt(this_data)
    station_flags = np.where(
        detrended_data > medfilt_nsig, True, False
    )

    # update flags with watershed algorithm
    station_flags = xrfi._ws_flag_waterfall(
        detrended_data, station_flags, nsig=ws_nsig
    )
    rfi_flags[this_slice] = station_flags
    detrended_rfi_data[this_slice] = detrended_data

# TODO: update the clobber part to be reasonable
clobber = True
history = "\nFlagged narrowband RFI using median filter detrending and "
history += f"a watershed algorithm. Median filter cut is {medfilt_nsig} " 
history += f"sigma; watershed algorithm cut is {ws_nsig}."
uvd.history += history

# Just write the detrended data to disk with the flags.
save_filename = os.path.join(
    save_dir, f"zen.{jd}.detrended.flagged.uvh5"
)
uvd.flag_array = rfi_flags
uvd.data_array = detrended_rfi_data
uvd.write_uvh5(save_filename, clobber=clobber)

