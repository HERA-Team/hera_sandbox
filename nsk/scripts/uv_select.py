#!/usr/bin/env python2.7
"""
read in miriad file, select on data and write out
"""
import pyuvdata
import numpy
import os
import glob
import argparse
import hera_cal as hc
from collections import OrderedDict as odict

args = argparse.ArgumentParser(description="read-in miriad file, select subset of data and write out")
args.add_argument("files", nargs='*', type=str, help="files to select data on")
args.add_argunent("--ext", default="S", type=str, help="file extension for output miriad file")
args.add_argument("--overwrite", default=False, action='store_true', help='overwrite output data files')
args.add_argument("--outdir", default=None, type=str, help="output directory")
args.add_argument("--ants", default=None, nargs='*', type=int, help="antenna numbers to keep")
args.add_argument("--freq_thin", default=None, type=int, help="factor by which to thin frequency axis")
args.add_argument("--time_thin", default=None, type=int, help="factor by which to thin time axis")
args.add_argument("--bl_types", default=None, nargs='*', type=tuple, help="(ant1, ant2) tuples specifying unique baseline types to keep")

# parse args
a = args.parse_args()

# iterate over files
for i, f in enumerate(a.files):
    # get outdir
    if a.outdir is None:
        outdir = os.path.dirname(f)
    else:
        outdir = a.outdir

    # check output
    fout = os.path.join(outdir, f + a.ext)
    if os.path.exists(fout) and a.overwrite is False:
        print "{} exists, not overwriting".format(fout)
        continue

    # Load into UVData object
    uvd = UVData()
    uvd.read_miriad(f)

    # get meta data
    freqs = uvd.freq_array.ravel()
    times = np.unique(uvd.time_array)
    antpos, ants = uvd.get_ENU_antpos(pick_data_ants=True)
    antpos = dict(zip(ants, antpos))
    bls = odict(map(lambda i: ((i[0], i[1]), None), zip(uvd.ant_1_array, uvd.ant_2_array))).keys()

    # select antennas
    if a.ants is not None:
        uvd.select(antenna_nums=a.ants)

    # select bls
    if a.bl_types is not None:
        reds = map(lambda blg: map(lambda bl: bl[:2], blg), hc.redcal.get_reds(antpos))
        keep = []
        for bl_group in reds:
            for bl in a.bl_types:
                if (bl in bl_group and bl in bls) or (bl[::-1] in bl_group and bl[::-1] in bls):
                    for b in bl_group:
                        if b in bls or b[::-1] in bls:
                            keep.append(b)

        uvd.select(ant_pairs_nums=keep)

    # select freqs
    if a.freq_thin is not None:
        uvd.select(frequencies=freqs[::a.freq_thin])
        uvd.channel_width = np.median(np.diff(np.unique(uvd.freq_array)))

    # select times
    if a.time_thin is not None:
        uvd.select(times=times[::a.time_thin])
        uvd.integration_time = np.median(np.diff(np.unique(uvd.time_array))) * 24 * 3600

    # write
    print "saving {}".format(fout)
    uvd.write_miriad(fout)





