#!/usr/bin/env python2.7
"""
sky_image.py
-------------
Visibility imaging with CASA 5.1.1

run this script as:
casa -c sky_image.py <args>

Nick Kern
nkern@berkeley.edu
Sept. 2018
"""
import sys
import os
import numpy as np
import argparse
import subprocess
import shutil
import glob

## Set Arguments
# Required Arguments
a = argparse.ArgumentParser(description="Run with casa as: casa -c sky_image.py <args>")
a.add_argument('--script', '-c', type=str, help='name of this script', required=True)
a.add_argument('--msin', default=None, type=str, help='path to a CASA measurement set. if fed a .uvfits, will convert to ms', required=True)
# IO Arguments
a.add_argument("--cleanspace", default=True, type=bool, help="Clean directory of image stem namespace before proceeding.")
a.add_argument('--source', default=None, type=str, help='Name of the main source in the field.')
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument("--silence", default=False, action='store_true', help="turn off output to stdout")
a.add_argument('--source_ext', default=None, type=str, help="Extension to default source name in output image files")
a.add_argument("--im_stem", default=None, type=str, help="Imagename stem for output images. Default is basename of input MS.")
# Imaging Arguments
a.add_argument('--image_mfs', default=False, action='store_true', help="make an MFS image across the selected band")
a.add_argument('--niter', default=[0], type=int, nargs='*', help='Number of clean iterations. Can be a list of niter for each mask provided.')
a.add_argument('--pxsize', default=300, type=int, help='pixel (cell) scale in arcseconds')
a.add_argument('--imsize', default=500, type=int, help='number of pixels along a side of the output square image.')
a.add_argument('--spw', default="", type=str, help="Imaging spectral window selection.")
a.add_argument('--uvrange', default="", type=str, help="CASA uvrange string to set in imaging.")
a.add_argument('--timerange', default=[""], type=str, nargs='*', help="Imaging timerange(s)")
a.add_argument('--spec_cube', default=False, action='store_true', help="image spectral cube as well as MFS.")
a.add_argument('--spec_dchan', default=40, type=int, help="number of channel averaging for a single image in the spectral cube.")
a.add_argument('--spec_start', default=100, type=int, help='starting channel for spectral cube')
a.add_argument('--spec_end', default=924, type=int, help='ending channel for spectral cube')
a.add_argument("--stokes", default='I', type=str, help="Polarizations to image. Cannot mix Stokes and Dipole pols. Ex. 'IQUV' or 'XXYY'. Default is 'I'")
a.add_argument("--mask", default=[''], type=str, nargs='*', help="CASA region string (or lists of them) to use as mask in CLEANing. Ex: 'circle[[1h55m0s,-30d40m0s],10deg]'")
a.add_argument("--weighting", default='briggs', type=str, help="Visibility weighting when imaging.")
a.add_argument("--robust", default=0, type=float, help="Robust parameter when briggs weighting.")
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--flag_autos', default=True, type=bool, help="flag autocorrelations in data.")
a.add_argument("--export_fits", default=True, type=bool, help="Export all output CASA images to FITS.")
# Plotting Arguments
a.add_argument("--plot_uvdist", default=False, action='store_true', help='make a uvdist plot')

def echo(message, type=0):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print("\n" + message + "\n" + "-"*40)

if __name__ == "__main__":
    # parse args
    args = a.parse_args()

    # get vars
    if args.source_ext is None:
        args.source_ext = ''
    verbose = args.silence is False

    # get phase center
    if args.source is not None:
        ra, dec = np.loadtxt('{}.loc'.format(args.source), dtype=str)
        ra, dec = ra.split(':'), dec.split(':')
        fixdir = 'J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec))
    else:
        args.source = ''
        fixdir = None

    msin = args.msin

    # get paths
    base_ms = os.path.basename(msin)
    if args.out_dir is None:
        out_dir = os.path.dirname(msin)
    else:
        out_dir = args.out_dir

    # check for uvfits
    if base_ms.split('.')[-1] == 'uvfits':
        echo("...converting uvfits to ms", type=1)
        uvfits = msin
        msin = os.path.join(out_dir, '.'.join(base_ms.split('.')[:-1] + ['ms']))
        base_ms = os.path.basename(msin)
        msfiles = glob.glob("{}*".format(msin))
        if len(msfiles) != 0:
            for i, msf in enumerate(msfiles):
                try:
                    os.remove(msf)
                except OSError:
                    shutil.rmtree(msf)
        importuvfits(uvfits, msin)
        echo("{}".format(msin))

    # rephase to source
    if fixdir is not None:
        echo("...fix vis to {} at {}".format(args.source, fixdir), type=1)
        fixvis(msin, msin, phasecenter=fixdir)

    # unflag
    if args.unflag is True:
        echo("...unflagging", type=1)
        flagdata(msin, mode='unflag')

    # flag autocorrs
    if args.flag_autos:
        echo("...flagging autocorrs", type=1)
        flagdata(msin, autocorr=True)

    # flag bad ants
    if args.ex_ants is not None:
        args.ex_ants = ','.join(map(lambda x: "HH"+x, args.ex_ants.split(',')))
        echo("...flagging bad ants: {}".format(args.ex_ants), type=1)
        flagdata(msin, mode='manual', antenna=args.ex_ants)

    # rflag
    if args.rflag is True:
        echo("...rfi flagging", type=1)
        flagdata(msin, mode='rflag')

    # get image stem
    if args.im_stem is None:
        im_stem = os.path.join(out_dir, base_ms + '.' + args.source + args.source_ext)
    else:
        im_stem = args.im_stem
  
    if args.cleanspace:
        # remove paths
        source_files = glob.glob(im_stem+'*')
        if len(source_files) > 0:
            for sf in source_files:
                if os.path.exists(sf):
                    try:
                        shutil.rmtree(sf)
                    except OSError:
                        os.remove(sf)

    # create mfs image
    if args.image_mfs:
        echo("...running MFS", type=1)
        def image_mfs(msin, im_stem, timerange, spw, niters, masks):
            assert len(niters) == len(masks), "len(niter) must equal len(mask)"
            for n, m in zip(niters, masks):
                echo("...cleaning {} for {} iters with mask '{}'".format(msin, n, m))
                clean(vis=msin, imagename=im_stem, spw=spw, niter=n, weighting=args.weighting, robust=args.robust, imsize=args.imsize,
                      cell='{}arcsec'.format(args.pxsize), mode='mfs', timerange=timerange, uvrange=args.uvrange, stokes=args.stokes,
                      mask=m)
            echo("...saving {}".format('{}.image'.format(im_stem)))
            if args.export_fits:
                exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))
                echo("...saving {}".format('{}.fits'.format(im_stem)))

        for i, tr in enumerate(args.timerange):
            if i == 0:
                image_mfs(msin, im_stem, tr, args.spw, args.niter, args.mask)
            else:
                image_mfs(msin, im_stem+'_tr{}'.format(i), tr, args.spw, args.niter, args.mask)

    # create spectrum
    if args.spec_cube:
        echo("...running MFS spectral cube clean", type=1)
        def spec_cube(msin, im_stem, timerange, spw, niters, masks):
            assert len(niters) == len(masks), "len(niter) must equal len(mask)"
            dchan = args.spec_dchan
            for i, chan in enumerate(np.arange(args.spec_start, args.spec_end, dchan)):
                spec_im_stem = '{}.spec{:04d}'.format(im_stem, chan)
                for n, m in zip(niters, masks):
                    echo("...cleaning {} for {} iters with mask '{}'".format(msin, n, m))
                    clean(vis=msin, imagename=spec_im_stem, niter=n, spw="0:{}~{}".format(chan, chan+dchan-1),
                            weighting=args.weighting, robust=args.robust, imsize=args.imsize, timerange=timerange, uvrange=args.uvrange,
                            cell='{}arcsec'.format(args.pxsize), mode='mfs', stokes=args.stokes, mask=m)
                echo("...saving {}.image".format(spec_im_stem))
                if args.export_fits:
                    exportfits(imagename='{}.image'.format(spec_im_stem), fitsimage='{}.fits'.format(spec_im_stem))
                    echo("...saving {}.fits".format(spec_im_stem))

        for i, tr in enumerate(args.timerange):
            if i == 0:
                spec_cube(msin, im_stem, tr, args.spw, args.niter, args.mask)
            else:
                spec_cube(msin, im_stem+'_tr{}'.format(i), tr, args.spw, args.niter, args.mask)

    # make uvdist plot
    if args.plot_uvdist:
        echo("...plotting uvdistance", type=1)
        # load visibility amplitudes
        ms.open(msin)
        data = ms.getdata(["amplitude", "antenna1", "antenna2", "uvdist", "axis_info", "flag"], ifraxis=True)
        amps = data['amplitude']
        flags = data['flag']
        uvdist = data['uvdist']
        freqs = data['axis_info']['freq_axis']['chan_freq']
        ms.close()
        # get rid of autos
        select = []
        for i, a1 in enumerate(data['antenna1']):
            if a1 != data['antenna2'][i]:
                select.append(i)
        amps = amps[:, :, select, :].squeeze()
        uvdist = uvdist[select, :]
        flags = flags[:, :, select, :].squeeze()
        # omit flagged data
        amps[flags] *= np.nan
        # average across time
        amps = np.nanmean(amps, axis=2)
        uvdist = np.nanmean(uvdist, axis=1)
        # average into channel bins
        freq_bins = np.median(np.split(np.linspace(100., 200., 1024, endpoint=True)[100:924], 4), axis=1)
        amps = np.nanmedian(np.split(amps[100:924, :], 4), axis=1)
        # plot
        import matplotlib.pyplot as plt
        def plot_uvdist(amp, ext):
            fig, ax = plt.subplots(1, 1, figsize=(10,7))
            ax.grid(True)
            p = ax.plot(uvdist, amp.T, marker='o', markersize=6, alpha=0.8, ls='')
            ax.set_xlabel("baseline length (meters)", fontsize=16)
            ax.set_ylabel("amplitude (arbitrary units)", fontsize=16)
            ax.set_title("{} for {}".format(ext, base_ms))
            ax.legend(p, map(lambda x: '{:0.0f} MHz'.format(x), freq_bins))
            file_fname = os.path.join(out_dir, "{}.{}.png".format(base_ms, ext))
            echo("...saving {}".format(file_fname))
            fig.savefig(file_fname, bbox_inches='tight', pad=0.05)
            plt.close()

        plot_uvdist(amps, "uvdist")


