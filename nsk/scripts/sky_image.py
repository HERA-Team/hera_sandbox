"""
sky_image.py
-------------
sky-based calibration with CASA 5.1.1

runt the script as:
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

# get args
a = argparse.ArgumentParser(description="Run with casa as: casa -c sky_image.py <args>")
a.add_argument('--script', '-c', type=str, help='name of this script', required=True)
a.add_argument('--msin', default=None, type=str, nargs='*', help='path to CASA measurement set(s). if fed a .uvfits, will convert to ms', required=True)
a.add_argument('--source', default=None, type=str, help='source name', required=True)
a.add_argument('--timerange', default="", type=str, help="calibration and clean timerange")
a.add_argument('--uvrange', default="", type=str, help="uvrange in meters (baseline length) to use in calibration and imaging")
a.add_argument('--refant', default=None, type=str, help='reference antenna')
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument('--nocal', default=False, action='store_true', help='skip calibration and just make an image')
a.add_argument('--noKGcal', default=False, action='store_true', help='do not perform K (dly) & G (phs) calibration')
a.add_argument('--noAcal', default=False, action='store_true', help='do not perform G (amp) calibration')
a.add_argument('--noBPcal', default=False, action='store_true', help='do not perform BandPass calibration (phs & amp)')
a.add_argument('--source_ext', default=None, type=str, help="extension to source name in output image files")
a.add_argument('--image_model', default=False, action='store_true', help="image model datacolumn instead of data datacolumn")
a.add_argument('--image_mfs', default=False, action='store_true', help="make an MFS image across the band")
a.add_argument('--spec_cube', default=False, action='store_true', help="image spectral cube as well as MFS.")
a.add_argument('--spec_dchan', default=40, type=int, help="number of channel averaging for a single image in the spectral cube.")
a.add_argument('--niter', default=50, type=int, help='number of clean iterations.')
a.add_argument('--pxsize', default=300, type=int, help='pixel (cell) scale in arcseconds')
a.add_argument('--imsize', default=500, type=int, help='number of pixels along a side of the output square image.')
a.add_argument('--bpoly', default=False, action='store_true', help="use BPOLY mode in bandpass")
a.add_argument('--degamp', default=4, type=int, help="amplitude polynomial degree for BPOLY")
a.add_argument('--degphase', default=1, type=int, help="phase polynomial degree for BPOLY")
a.add_argument('--cleanspw', default="0:200~850", type=str, help="spectral window selection for clean")
a.add_argument("--plot_uvdist", default=False, action='store_true', help='make a uvdist plot')
a.add_argument("--silence", default=False, action='store_true', help="turn off output to stdout")

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
    KGcal = args.noKGcal is False
    Acal = args.noAcal is False
    BPcal = args.noBPcal is False

    # check for .loc and .cl files
    if os.path.exists('{}.loc'.format(args.source)) is False:
        raise AttributeError("{}.loc file doesn't exist in working directory".format(args.source))
    if os.path.exists('{}.cl'.format(args.source)) is False:
        raise AttributeError("{}.cl file doesn't exist in working directory".format(args.source))

    # configure refant
    if args.refant is None and args.nocal is False:
        raise AttributeError("if nocal is False, refant needs to be specified")
    args.refant = "HH" + str(args.refant)

    # get phase center
    ra, dec = np.loadtxt('{}.loc'.format(args.source), dtype=str)
    ra, dec = ra.split(':'), dec.split(':')
    fixdir = 'J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec))
    msin = args.msin
    if len(msin) == 1:
        # get paths
        msin = msin[0]
    else:
        # iterate through files and convert if uvfits
        ms_list = []
        for i, m in enumerate(msin):
            if m.split('.')[-1] == 'uvfits':
                uvfms = '.'.join(m.split('.')[:-1] + ["ms"])
                if os.path.exists(uvfms):
                    shutil.rmtree(uvfms)
                importuvfits(m, uvfms)
                ms_list.append(uvfms)
            elif m.split('.')[-1] == 'ms':
                ms_list.append(m)
        # basename and path is first ms in ms_list
        msin = ms_list[0]

        # rephase to source
        echo("...fix vis to {}".format(fixdir), type=1)
        for m in ms_list:
            fixvis(m, m, phasecenter=fixdir)

        # concatenate all ms into zeroth ms
        echo("...concatenating visibilities", type=1)
        concat(vis=ms_list, concatvis=msin, timesort=True)

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

    # rephase to source
    echo("...fix vis to {}".format(fixdir), type=1)
    fixvis(msin, msin, phasecenter=fixdir)

    # insert source model
    if args.nocal is False or args.image_model is True:
        ft(msin, complist="{}.cl".format(args.source), usescratch=True)

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
        gaincal(msin, caltable=kc, gaintype="K", solint='inf', refant=args.refant, minsnr=3.0, spw="0:100~924",
                gaintable=gaintables, timerange=args.timerange, uvrange=args.uvrange)
        plotcal(kc, xaxis='antenna', yaxis='delay', figfile='{}.png'.format(kc), showgui=False)
        gaintables.append(kc)

        # write delays to text file
        tb.open(kc)
        delays = tb.getcol('FPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]]
        delay_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
        tb.close()
        np.savetxt("{}.csv".format(kc), np.vstack([delay_ants, delays]).T, fmt="%12.8f", header="ant_num, delay [ns]", delimiter=", ")
        echo("...Saving delays to {}.csv".format(kc))
        echo("...Saving plotcal to {}.png".format(kc))

        # perform initial G calibration for phase (per-spw and per-pol gain)
        echo("...performing G gaincal for phase", type=1)
        if os.path.exists(gpc):
            shutil.rmtree(gpc)
        if os.path.exists("{}.png".format(gpc)):
            os.remove("{}.png".format(gpc))
        gaincal(msin, caltable=gpc, gaintype='G', solint='inf', refant=args.refant, minsnr=3, calmode='p',
                spw="0:100~924", gaintable=gaintables, timerange=args.timerange, uvrange=args.uvrange)
        plotcal(gpc, xaxis='antenna', yaxis='phase', figfile='{}.png'.format(gpc), showgui=False)
        gaintables.append(gpc)

        # write phase to file
        tb.open(gpc)
        phases = np.angle(tb.getcol('CPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]])
        phase_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
        tb.close()
        np.savetxt("{}.csv".format(gpc), np.vstack([phase_ants, phases]).T, fmt="%12.8f", header="ant_num, phase [radians]", delimiter=", ")
        echo("...Saving phases to {}.csv".format(gpc))
        echo("...Saving plotcal to {}.png".format(gpc))

        return gaintables

    def ACAL(msin, gaintables=[]):
        # gaincal G amplitude
        echo("...performing G gaincal for amplitude", type=1)
        gac = msin+'.{}.cal'.format('G')
        if os.path.exists(gac):
            shutil.rmtree(gac)
        if os.path.exists("{}.png".format(gac)):
            os.remove("{}.png".format(gac))
        gaincal(msin, caltable=gac, gaintype='G', solint='inf', refant=args.refant, minsnr=3, calmode='a',
                spw="0:100~924", gaintable=gaintables, timerange=args.timerange, uvrange=args.uvrange)
        plotcal(gac, xaxis='antenna', yaxis='amp', figfile='{}.png'.format(gac), showgui=False)
        gaintables.append(gac)

        # write amp to file
        tb.open(gac)
        amps = np.abs(tb.getcol('CPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]])
        amp_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
        tb.close()
        np.savetxt("{}.csv".format(gac), np.vstack([amp_ants, amps]).T, fmt="%12.8f", header="ant_num, amplitude", delimiter=", ")
        echo("...Saving amps to {}.csv".format(gac))
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
        bandpass(vis=msin, spw="0:100~924", minsnr=2, bandtype=Btype, degamp=args.degamp, degphase=args.degphase,
                caltable=bc, gaintable=gaintables, solint='inf', refant=args.refant, timerange=args.timerange, uvrange=args.uvrange)
        plotcal(bc, xaxis='chan', yaxis='amp', figfile="{}.amp.png".format(bc), showgui=False)
        plotcal(bc, xaxis='chan', yaxis='phase', figfile="{}.phs.png".format(bc), showgui=False)
        gaintables.append(bc)

        # write bp to file
        if args.bpoly is False:
            # get flags and bandpass data
            tb.open(bc)
            flags = ~tb.getcol('FLAG')[0]
            flagged_ants = np.sum(flags, axis=0).astype(np.bool)
            flags = flags[:, flagged_ants]
            bp = tb.getcol('CPARAM')[0][:, flagged_ants]
            bp_ants = tb.getcol("ANTENNA1")[flagged_ants]
            tb.close()
            # load spectral window data
            tb.open(bc+"/SPECTRAL_WINDOW")
            freqs = tb.getcol("CHAN_FREQ")
            tb.close()
            # write to file
            bp_data = np.concatenate([freqs/1e6, bp.real, bp.imag, ~flags], axis=1)
            bp_ants_str = ', '.join(map(lambda x: str(x)+'r', bp_ants)) + \
                            ', ' + ', '.join(map(lambda x: str(x)+'i', bp_ants)) + \
                            ', ' + ', '.join(map(lambda x: str(x)+'f', bp_ants))
            np.savetxt("{}.csv".format(bc), bp_data, fmt="%12.8f", header="freq (MHz), {}".format(bp_ants_str), delimiter=", ")
            echo("...Saving bandpass to {}.csv".format(bc))
            echo("...Saving amp plotcal to {}.amp.png".format(bc))
            echo("...Saving phs plotcal to {}.phs.png".format(bc))
        else:
            echo("couldn't access bandpass data. Note BPOLY solutions not currently compatible w/ caltable2calfits.py")

        return gaintables

    if args.nocal is False:
        ## Begin Calibration ##
        # run through various calibration options
        gaintables = []
        if KGcal:
            gaintables = KGCAL(msin, gaintables)

        if Acal:
            gaintables = ACAL(msin, gaintables)

        if BPcal:
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
        ms_split = msin

    if args.image_mfs is True or args.spec_cube is True:
        # clean
        im_stem = os.path.join(out_dir, base_ms + '.' + args.source + args.source_ext)
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
        model_im_stem = os.path.join(out_dir, base_ms + '.model.' + args.source + args.source_ext)
        split(ms_split, model_ms_split, datacolumn='model')


    # create mfs image
    if args.image_mfs:
        echo("...running MFS clean", type=1)
        if args.cleanspw is None:
            cleanspw = "0:100~924"
        else:
            cleanspw = args.cleanspw

        def image_mfs(msin, im_stem):
            clean(vis=msin, imagename=im_stem, spw=cleanspw, niter=args.niter, weighting='briggs', robust=0, imsize=[args.imsize, args.imsize],
                  cell=['{}arcsec'.format(args.pxsize)], mode='mfs', timerange=args.timerange, uvrange=args.uvrange)
            exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))
            print("...saving {}".format('{}.fits'.format(im_stem)))

        image_mfs(ms_split, im_stem)
   
        if args.image_model:
            image_mfs(model_ms_split, model_im_stem)

    # create spectrum
    if args.spec_cube:
        echo("...running spectral cube clean", type=1)
        def spec_cube(msin, im_stem):
            dchan = args.spec_dchan
            for i, chan in enumerate(np.arange(100, 924, dchan)):
                clean(vis=msin, imagename=im_stem+'.spec{}'.format(chan), niter=0, spw="0:{}~{}".format(chan, chan+dchan-1),
                        weighting='briggs', robust=0, imsize=[args.imsize, args.imsize], timerange=args.timerange, uvrange=args.uvrange,
                        cell=['{}arcsec'.format(args.pxsize)], mode='mfs')#, mask='circle[[{}h{}m{}s, {}d{}m{}s ], 7deg]'.format(*(ra+dec)))
                exportfits(imagename='{}.spec{}.image'.format(im_stem, chan), fitsimage='{}.spec{}.fits'.format(im_stem, chan))
                print("...saving {}".format('{}.spec{}.fits'.format(im_stem, chan)))

        spec_cube(ms_split, im_stem)

        if args.image_model:
            spec_cube(model_ms_split, model_im_stem)

    # make uvdist plot
    if args.plot_uvdist:
        echo("...plotting uvdistance")
        # add model if nocal
        if args.nocal:
            ft(ms_split, complist="{}.cl".format(args.source), usescratch=True)
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
        amps[~flags] *= np.nan
        mamps[~flags] *= np.nan
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



