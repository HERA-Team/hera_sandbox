#!/usr/bin/env python2.7
"""
sky_image.py
-------------
sky-based calibration with CASA 5.1.1

run the script as:
casa -c sky_image.py <args>

Nick Kern
nkern@berkeley.edu
Nov. 2017
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
a.add_argument('--source', default=None, type=str, help='source name', required=True)
# IO Arguments
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument("--silence", default=False, action='store_true', help="turn off output to stdout")
a.add_argument('--source_ext', default=None, type=str, help="Extension to default source name in output image files")
a.add_argument("--im_stem", default=None, type=str, help="Imagename stem for output images. Default is basename of input MS.")
# Calibration Arguments
a.add_argument("--model_im", default=None, type=str, help="path to model image, if None will look for a {source}.cl file")
a.add_argument('--refant', default=None, type=str, help='reference antenna')
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--KGcal', default=False, action='store_true', help='perform K (dly) & G (phs) calibration')
a.add_argument("--KGsnr", default=2.0, type=float, help="KG calibration Signal-to-Noise cut")
a.add_argument('--Acal', default=False, action='store_true', help='perform G (amp) calibration')
a.add_argument("--Asnr", default=2.0, type=float, help="G-amplitude calibration Signal-to-Noise cut")
a.add_argument('--BPcal', default=False, action='store_true', help='perform BandPass calibration (phs & amp)')
a.add_argument("--BPsnr", default=2.0, type=float, help="bandpass calibration Signal-to-Noise cut")
a.add_argument('--uvrange', default="", type=str, help="uvrange in meters (baseline length) to use in calibration and imaging")
a.add_argument('--timerange', default=[""], type=str, nargs='*', help="calibration and clean timerange(s)")
a.add_argument('--bpoly', default=False, action='store_true', help="use BPOLY mode in bandpass")
a.add_argument('--degamp', default=4, type=int, help="amplitude polynomial degree for BPOLY")
a.add_argument('--degphase', default=1, type=int, help="phase polynomial degree for BPOLY")
a.add_argument('--calspw', default='0:100~924', type=str, help="Calibration spectral window selection")
a.add_argument('--smodel', default=[], type=float, nargs='*', help="Stokes source model as I Q U V")
# Imaging Arguments
a.add_argument('--image_mfs', default=False, action='store_true', help="make an MFS image across the band")
a.add_argument('--niter', default=50, type=int, help='number of clean iterations.')
a.add_argument('--pxsize', default=300, type=int, help='pixel (cell) scale in arcseconds')
a.add_argument('--imsize', default=500, type=int, help='number of pixels along a side of the output square image.')
a.add_argument('--cleanspw', default="0:100~924", type=str, help="spectral window selection for clean")
a.add_argument('--image_model', default=False, action='store_true', help="image model datacolumn instead of data datacolumn")
a.add_argument('--spec_cube', default=False, action='store_true', help="image spectral cube as well as MFS.")
a.add_argument('--spec_dchan', default=40, type=int, help="number of channel averaging for a single image in the spectral cube.")
a.add_argument('--spec_start', default=100, type=int, help='starting channel for spectral cube')
a.add_argument('--spec_end', default=924, type=int, help='ending channel for spectral cube')
a.add_argument("--stokes", default='I', type=str, help="Stokes parameters to image.")
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

    # check for .loc and .cl files
    if os.path.exists('{}.loc'.format(args.source)) is False:
        raise AttributeError("{}.loc file doesn't exist in working directory".format(args.source))

    # configure refant
    if args.refant is None and (args.KGcal is True or args.Acal is True or args.BPcal is True):
        raise AttributeError("if calibrating, refant needs to be specified")
    args.refant = "HH" + str(args.refant)

    # get phase center
    ra, dec = np.loadtxt('{}.loc'.format(args.source), dtype=str)
    ra, dec = ra.split(':'), dec.split(':')
    fixdir = 'J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec))
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
    echo("...fix vis to {}".format(fixdir), type=1)
    fixvis(msin, msin, phasecenter=fixdir)

    # insert source model
    if (args.KGcal is True or args.Acal is True or args.BPcal is True) or args.image_model is True:
        if args.model_im is None:
            echo("...inserting {} as MODEL".format("{}.cl".format(args.source)), type=1)
            ft(msin, complist="{}.cl".format(args.source), usescratch=True)
        else:
            echo("...inserting {} as MODEL".format(args.model_im), type=1)
            ft(msin, model=args.model_im, usescratch=True)

    # unflag
    if args.unflag is True:
        echo("...unflagging", type=1)
        flagdata(msin, mode='unflag')

    # flag autocorrs
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

    def KGCAL(msin, gaintables=[]):
        ## perform per-antenna delay and phase calibration ##
        # setup calibration tables     
        kc = os.path.join(out_dir, base_ms+'.{}.cal'.format('K'))
        gpc = os.path.join(out_dir, base_ms+'.{}.cal'.format('Gphs'))

        # perform initial K calibration (per-antenna delay)
        echo("...performing K gaincal", type=1)
        if os.path.exists(kc):
            shutil.rmtree(kc)
        if os.path.exists("{}.png".format(kc)):
            os.remove("{}.png".format(kc))
        gaincal(msin, caltable=kc, gaintype="K", solint='inf', refant=args.refant, minsnr=args.KGsnr, spw=args.calspw,
                gaintable=gaintables, timerange=cal_timerange, uvrange=args.uvrange)
        plotcal(kc, xaxis='antenna', yaxis='delay', figfile='{}.png'.format(kc), showgui=False)
        gaintables.append(kc)

        # write delays as npz file
        tb.open(kc)
        delays = tb.getcol('FPARAM')[0, 0]
        delay_ants = tb.getcol('ANTENNA1')
        delay_flags = tb.getcol('FLAG')[0, 0]
        tb.close()
        np.savez("{}.npz".format(kc), delay_ants=delay_ants, delays=delays, delay_flags=delay_flags)
        echo("...Saving delays to {}.npz".format(kc))
        echo("...Saving plotcal to {}.png".format(kc))

        # perform initial G calibration for phase (per-spw and per-pol gain)
        echo("...performing G gaincal for phase", type=1)
        if os.path.exists(gpc):
            shutil.rmtree(gpc)
        if os.path.exists("{}.png".format(gpc)):
            os.remove("{}.png".format(gpc))
        gaincal(msin, caltable=gpc, gaintype='G', solint='inf', refant=args.refant, minsnr=args.KGsnr, calmode='p',
                spw=args.calspw, gaintable=gaintables, timerange=cal_timerange, uvrange=args.uvrange)
        plotcal(gpc, xaxis='antenna', yaxis='phase', figfile='{}.png'.format(gpc), showgui=False)
        gaintables.append(gpc)

        # write phase to file
        tb.open(gpc)
        phases = np.angle(tb.getcol('CPARAM')[0, 0])
        phase_ants = tb.getcol('ANTENNA1')
        phase_flags = tb.getcol('FLAG')[0, 0]
        tb.close()
        np.savez("{}.npz".format(gpc), phase_ants=phase_ants, phases=phases, phase_flags=phase_flags)
        echo("...Saving phases to {}.npz".format(gpc))
        echo("...Saving plotcal to {}.png".format(gpc))

        return gaintables

    def ACAL(msin, gaintables=[]):
        # gaincal G amplitude
        echo("...performing G gaincal for amplitude", type=1)
        gac = msin+'.{}.cal'.format('Gamp')
        if os.path.exists(gac):
            shutil.rmtree(gac)
        if os.path.exists("{}.png".format(gac)):
            os.remove("{}.png".format(gac))
        gaincal(msin, caltable=gac, gaintype='G', solint='inf', refant=args.refant, minsnr=args.Asnr, calmode='a',
                spw=args.calspw, gaintable=gaintables, timerange=cal_timerange, uvrange=args.uvrange)
        plotcal(gac, xaxis='antenna', yaxis='amp', figfile='{}.png'.format(gac), showgui=False)
        gaintables.append(gac)

        # write amp to file
        tb.open(gac)
        amps = np.abs(tb.getcol('CPARAM')[0, 0])
        amp_ants = tb.getcol('ANTENNA1')
        amp_flags = tb.getcol('FLAG')[0, 0]
        tb.close()
        np.savez("{}.npz".format(gac), amp_ants=amp_ants, amps=amps, amp_flags=amp_flags)
        echo("...Saving amps to {}.npz".format(gac))
        echo('...Saving G amp plotcal to {}.png'.format(gac))

        return gaintables

    def BPCAL(msin, gaintables=[]):
        # calibrated bandpass
        echo("...performing B bandpass cal", type=1)
        bc = msin+'.{}.cal'.format('B')
        Btype = "B"
        if args.bpoly:
            Btype="BPOLY"
        if os.path.exists(bc):
            shutil.rmtree(bc)
        if os.path.exists("{}.amp.png".format(bc)):
            os.remove("{}.amp.png".format(bc))
        if os.path.exists("{}.phs.png".format(bc)):
            os.remove("{}.phs.png".format(bc))
        bandpass(vis=msin, spw="", minsnr=args.BPsnr, bandtype=Btype, degamp=args.degamp, degphase=args.degphase,
                caltable=bc, gaintable=gaintables, solint='inf', refant=args.refant, timerange=cal_timerange,
                uvrange=args.uvrange, smodel=args.smodel)
        plotcal(bc, xaxis='chan', yaxis='amp', figfile="{}.amp.png".format(bc), showgui=False)
        plotcal(bc, xaxis='chan', yaxis='phase', figfile="{}.phs.png".format(bc), showgui=False)
        gaintables.append(bc)

        # write bp to file
        if args.bpoly is False:
            # get flags and bandpass data
            tb.open(bc)
            bp = tb.getcol('CPARAM')[0]
            bp_ants = tb.getcol("ANTENNA1")
            bp_flags = tb.getcol('FLAG')[0]
            tb.close()
            # load spectral window data
            tb.open(bc+"/SPECTRAL_WINDOW")
            bp_freqs = tb.getcol("CHAN_FREQ")
            tb.close()
            # write to file
            np.savez("{}.npz".format(bc), bp=bp, bp_ants=bp_ants, bp_flags=bp_flags, bp_freqs=bp_freqs)
            echo("...Saving bandpass to {}.npz".format(bc))
            echo("...Saving amp plotcal to {}.amp.png".format(bc))
            echo("...Saving phs plotcal to {}.phs.png".format(bc))
        else:
            echo("couldn't access bandpass data. Note BPOLY solutions not currently compatible w/ caltable2calfits.py")

        return gaintables

    if (args.KGcal is True or args.Acal is True or args.BPcal is True):
        ## Begin Calibration ##
        # init cal_timerange
        cal_timerange = ','.join(args.timerange)
        # run through various calibration options
        gaintables = []
        if args.KGcal:
            gaintables = KGCAL(msin, gaintables)

        if args.Acal:
            gaintables = ACAL(msin, gaintables)

        if args.BPcal:
            gaintables = BPCAL(msin, gaintables)

        # apply calibration gaintables
        echo("...applying gaintables: \n {}".format('\n'.join(gaintables)), type=1)
        applycal(msin, gaintable=gaintables)
        ms_split = os.path.join(out_dir, "{}.split".format(base_ms))
        files = glob.glob("{}*".format(ms_split))
        for f in files:
            if os.path.exists(f):
                try:
                    shutil.rmtree(f)
                except OSError:
                    os.remove(f)

        split(msin, ms_split, datacolumn="corrected")
        fixvis(ms_split, ms_split, phasecenter=fixdir)

    else:
        echo("...no calibration performed", type=1)
        ms_split = msin

    if args.image_mfs is True or args.spec_cube is True:
        # remove paths
        if args.im_stem is None:
            im_stem = os.path.join(out_dir, base_ms + '.' + args.source + args.source_ext)
        else:
            im_stem = args.im_stem
        echo("...performing clean for output files:\n{}.*".format(im_stem), type=1)
        source_files = glob.glob(im_stem+'*')
        if len(source_files) > 0:
            for sf in source_files:
                if os.path.exists(sf):
                    try:
                        shutil.rmtree(sf)
                    except OSError:
                        os.remove(sf)

    if args.image_model:
        model_ms_split = ms_split + ".model"
        if args.im_stem is None:
            model_im_stem = os.path.join(out_dir, base_ms + '.model.' + args.source + args.source_ext)
        else:
            model_im_stem = args.im_stem + '.model'
        if args.model_im is None:
            echo("...inserting {} as MODEL".format("{}.cl".format(args.source)), type=1)
            ft(ms_split, complist="{}.cl".format(args.source), usescratch=True)
        else:
            echo("...inserting {} as MODEL".format(args.model_im), type=1)
            ft(ms_split, model=args.model_im, usescratch=True)
        split(ms_split, model_ms_split, datacolumn='model')

    # create mfs image
    if args.image_mfs:
        echo("...running MFS clean", type=1)
        def image_mfs(msin, im_stem, timerange, cleanspw):
            clean(vis=msin, imagename=im_stem, spw=cleanspw, niter=args.niter, weighting='briggs', robust=0, imsize=[args.imsize, args.imsize],
                  cell=['{}arcsec'.format(args.pxsize)], mode='mfs', timerange=timerange, uvrange=args.uvrange, stokes=args.stokes)
            exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))
            print("...saving {}".format('{}.fits'.format(im_stem)))

    
        for i, tr in enumerate(args.timerange):
            if i == 0:
                image_mfs(ms_split, im_stem, tr, args.cleanspw)
                if args.image_model:
                    image_mfs(model_ms_split, model_im_stem, tr, args.cleanspw)
            else:
                image_mfs(ms_split, im_stem+'_tr{}'.format(i), tr, args.cleanspw)
                if args.image_model:
                    image_mfs(model_ms_split, model_im_stem+'_tr{}'.format(i), tr, args.cleanspw)

    # create spectrum
    if args.spec_cube:
        echo("...running spectral cube clean", type=1)
        def spec_cube(msin, im_stem, timerange):
            dchan = args.spec_dchan
            for i, chan in enumerate(np.arange(args.spec_start, args.spec_end, dchan)):
                clean(vis=msin, imagename=im_stem+'.spec{:04d}'.format(chan), niter=0, spw="0:{}~{}".format(chan, chan+dchan-1),
                        weighting='briggs', robust=0, imsize=[args.imsize, args.imsize], timerange=timerange, uvrange=args.uvrange,
                        cell=['{}arcsec'.format(args.pxsize)], mode='mfs', stokes=args.stokes)#, mask='circle[[{}h{}m{}s, {}d{}m{}s ], 7deg]'.format(*(ra+dec)))
                exportfits(imagename='{}.spec{:04d}.image'.format(im_stem, chan), fitsimage='{}.spec{:04d}.fits'.format(im_stem, chan))
                print("...saving {}".format('{}.spec{:04d}.fits'.format(im_stem, chan)))

        for i, tr in enumerate(args.timerange):
            if i == 0:
                spec_cube(ms_split, im_stem, tr)
                if args.image_model:
                    spec_cube(model_ms_split, model_im_stem, tr)
            else:
                spec_cube(ms_split, im_stem+'_tr{}'.format(i), tr)
                if args.image_model:
                    spec_cube(model_ms_split, model_im_stem+'_tr{}'.format(i), tr)

    # make uvdist plot
    if args.plot_uvdist:
        echo("...plotting uvdistance", type=1)
        # add model to ms_split
        if args.model_im is None:
            ft(ms_split, complist="{}.cl".format(args.source), usescratch=True)
        else:
            ft(ms_split, model=args.model_im, usescratch=True)
        # load visibility amplitudes
        ms.open(ms_split)
        data = ms.getdata(["amplitude", "antenna1", "antenna2", "uvdist", "axis_info", "flag", "model_amplitude"], ifraxis=True)
        amps = data['amplitude']
        mamps = data['model_amplitude']
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
        mamps = mamps[:, :, select, :].squeeze()
        uvdist = uvdist[select, :]
        flags = flags[:, :, select, :].squeeze()
        # omit flagged data
        amps[flags] *= np.nan
        mamps[flags] *= np.nan
        # average across time
        amps = np.nanmean(amps, axis=2)
        mamps = np.nanmean(mamps, axis=2)
        uvdist = np.nanmean(uvdist, axis=1)
        # average into channel bins
        freq_bins = np.median(np.split(np.linspace(100., 200., 1024, endpoint=True)[100:924], 4), axis=1)
        amps = np.nanmedian(np.split(amps[100:924, :], 4), axis=1)
        mamps = np.nanmedian(np.split(mamps[100:924, :], 4), axis=1)
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

        plot_uvdist(amps, "data_uvdist")
        plot_uvdist(mamps, "model_uvdist")



