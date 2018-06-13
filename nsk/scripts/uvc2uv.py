#!/usr/bin/env python2.7
"""
uvc2uv.py
---------------

convert *.uvc files to *.uv files

Nicholas Kern
Jan. 2018
"""
from pyuvdata import UVData
import numpy as np
import argparse
import os
import copy
import aipy
import hera_cal as hc
from pyuvdata import uvutils

args = argparse.ArgumentParser(description="")

args.add_argument("files", type=str, nargs='*', help='path to *.uvc files to convert to *.uv files')
args.add_argument("--calfile", type=str, help="path to calfile stem found in user's PATH")
args.add_argument("--outdir", default=None, type=str, help="output directory")
args.add_argument("--overwrite", default=False, action='store_true', help="overwrite output files")
args.add_argument("--name_prefix", default="HH", help="prefix to antenna number for antenna name")

def uvc2uv(uvcfile, calfile, outdir=None, overwrite=False, name_prefix="HH"):
    """
    """
    # check output
    if outdir is None:
        outdir = os.path.dirname(uvcfile)

    uvfile = os.path.splitext(os.path.basename(uvcfile))
    uvfile = uvfile[0] + ".uv"
    output_fname = os.path.join(outdir, uvfile)
    if os.path.exists(output_fname) and overwrite is False:
        raise IOError("...{} exists, not overwriting")

    # load file
    uvc = UVData()
    uvc.read_miriad(uvcfile)

    # get antenna numbers from data
    uvc_ants = np.unique(np.concatenate([uvc.ant_1_array, uvc.ant_2_array]))

    # reorder according to data
    antenna_numbers = copy.copy(uvc.antenna_numbers)
    uvc_antcut = [True if a in uvc_ants else False for a in antenna_numbers]
    antenna_numbers = antenna_numbers[uvc_antcut]
    uvc_antsort = []
    for a in uvc_ants:
        if a in antenna_numbers:
            uvc_antsort.append(antenna_numbers.tolist().index(a))
    uvc_ants = uvc_ants[uvc_antsort]
    Nants_data = len(uvc_ants)
    if uvc.antenna_diameters is not None:
        uvc_ant_diameters = uvc.antenna_diameters[uvc_antcut][uvc_antsort]
    else:
        uvc_ant_diameters = np.ones_like(uvc_antcut, np.float)

    # get antenna positions from calfile
    aa = aipy.cal.get_aa(calfile, np.array([0.15]))
    info = hc.omni.aa_to_info(aa, pols=['y'], tol=5.0)
    cf_ants = info.subsetant
    cf_antpos = info.antloc

    # check ants in uvc are at least in cf
    assert len(set(uvc_ants) - set(cf_ants)) == 0, "ants {} found in data but not calfile".format(set(uvc_ants) - set(cf_ants))

    # reorder cf_ants to match uvc_ants
    cf_antsort = [cf_ants.tolist().index(a) for a in uvc_ants]
    cf_ants = cf_ants[cf_antsort]
    cf_antpos = cf_antpos[cf_antsort]

    # convert antpos from ENU meters to ITRF nanosec
    ECEF_antpos = uvutils.ECEF_from_ENU(cf_antpos.T, *uvc.telescope_location_lat_lon_alt).T - uvc.telescope_location

    # reset antenna position and number information
    uvc.Nants_data = Nants_data
    uvc.Nants_telescope = Nants_data
    uvc.antenna_numbers = cf_ants
    uvc.antenna_names = map(lambda a: "{}{}".format(name_prefix, a), cf_ants)
    uvc.antenna_positions = ECEF_antpos
    uvc.antenna_diameters = uvc_ant_diameters

    # write to file
    print("..saving {}".format(output_fname))
    uvc.write_miriad(output_fname, clobber=True)

if __name__ == "__main__":

    # parse args
    a = args.parse_args()
    kwargs = dict(vars(a))
    kwargs.pop('files')
    kwargs.pop('calfile')

    # iterate over files
    for i, f in enumerate(a.files):
        print("...loading {}".format(f))
        uvc2uv(f, a.calfile, **kwargs)

