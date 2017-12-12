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
a.add_argument('--refant', default=None, type=str, help='reference antenna')
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument('--nocal', default=False, action='store_true', help='skip calibration and just make an image')
a.add_argument('--onlyKcal', default=False, action='store_true', help='only perform K calibration, then image')
a.add_argument('--onlyKGcal', default=False, action='store_true', help='only perform K & G (phs) calibration, then image')
a.add_argument('--source_ext', default=None, type=str, help="extension to source name in output image files")
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

args = a.parse_args()

if __name__ == "__main__":
    # get vars
    if args.source_ext is None:
        source_ext = ''

    # configure refant
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
        print("...fix vis to {}".format(fixdir))
        for m in ms_list:
            fixvis(m, m, phasecenter=fixdir)

        # concatenate all ms into zeroth ms
        print("...concatenating visibilities")
        concat(vis=ms_list, concatvis=msin, timesort=True)

    # get paths
    base_ms = os.path.basename(msin)
    if args.out_dir is None:
        out_dir = os.path.dirname(msin)
    else:
        out_dir = args.out_dir

    # check for uvfits
    if base_ms.split('.')[-1] == 'uvfits':
        print("...converting uvfits to ms")
        uvfits = msin
        msin = os.path.join(out_dir, '.'.join(base_ms.split('.')[:-1] + ['ms']))
        base_ms = os.path.basename(msin)
        if os.path.exists(msin):
            shutil.rmtree(msin)
        importuvfits(uvfits, msin)

    # rephase to source
    print("...fix vis to {}".format(fixdir))
    fixvis(msin, msin, phasecenter=fixdir)

    # insert source model
    ft(msin, complist="{}.cl".format(args.source), usescratch=True)

    # unflag
    if args.unflag is True:
        print("...unflagging")
        flagdata(msin, mode='unflag')

    # flag autocorrs
    print("...flagging autocorrs")
    flagdata(msin, autocorr=True)

    # flag bad ants
    if args.ex_ants is not None:
        args.ex_ants = ','.join(map(lambda x: "HH"+x, args.ex_ants.split(',')))
        print("...flagging bad ants: {}".format(args.ex_ants))
        flagdata(msin, mode='manual', antenna=args.ex_ants)

    # rflag
    if args.rflag is True:
        print("...rfi flagging")
        flagdata(msin, mode='rflag')

    # setup calibration tables
    kc = os.path.join(out_dir, base_ms+'.{}.cal'.format('K'))
    gc = os.path.join(out_dir, base_ms+'.{}.cal'.format('G'))

    if args.nocal is True:
        Bsplit = msin
    else:
        # perform initial K calibration (per-antenna delay)
        print("...performing K gaincal")
        if os.path.exists(kc):
            shutil.rmtree(kc)
        if os.path.exists("{}.png".format(kc)):
            os.remove("{}.png".format(kc))
        gaincal(msin, caltable=kc, gaintype="K", solint='inf', refant=args.refant, minsnr=3.0, spw="0:100~924",
                gaintable=[], timerange=args.timerange)
        plotcal(kc, xaxis='antenna', yaxis='delay', figfile='{}.png'.format(kc), showgui=False)

        # write delays to text file
        tb.open(kc)
        delays = tb.getcol('FPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]]
        delay_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
        tb.close()
        np.savetxt("{}.csv".format(kc), np.vstack([delay_ants, delays]).T, fmt="%12.8f", header="ant_num, delay [ns]", delimiter=", ")
        print("...Saving delays to {}.csv".format(kc))

        if args.onlyKcal:
            applycal(msin, gaintable=[kc])

        else:
            # perform initial G calibration for phase (per-spw and per-pol gain)
            print("...performing G gaincal for phase")
            if os.path.exists(gc):
                shutil.rmtree(gc)
            gaincal(msin, caltable=gc, gaintype='G', solint='inf', refant=args.refant, minsnr=3, calmode='p',
                    spw="0:100~924", gaintable=[kc], timerange=args.timerange)
            plotcal(gc, xaxis='antenna', yaxis='phase', figfile='{}.png'.format(gc), showgui=False)

            # write phase to file
            tb.open(gc)
            phases = np.angle(tb.getcol('CPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]])
            phase_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
            tb.close()
            np.savetxt("{}.csv".format(gc), np.vstack([phase_ants, phases]).T, fmt="%12.8f", header="ant_num, phase [radians]", delimiter=", ")
            print("...Saving phases to {}.csv".format(gc))

            # apply calibrations
            print("...applying KG gaincal")
            applycal(msin, gaintable=[kc, gc])

        # split data
        KGsplit = os.path.join(out_dir, "{}.KGsplit".format(base_ms))
        files = glob.glob("{}*".format(KGsplit))
        for f in files:
            if os.path.exists(f):
                try:
                    shutil.rmtree(f)
                except OSError:
                    os.remove(f)
        split(msin, KGsplit, datacolumn="corrected")
        fixvis(KGsplit, KGsplit, phasecenter=fixdir)
        ft(KGsplit, complist="{}.cl".format(args.source), usescratch=True)

        if args.onlyKGcal is True or args.onlyKcal is True:
            Bsplit = KGsplit
        else:
            # calibrated bandpass
            print("...performing B bandpass cal")
            bc = KGsplit+'.{}.cal'.format('B')
            Btype = "B"
            if args.bpoly:
                Btype="BPOLY"
            bandpass(vis=KGsplit, spw="0:100~924", minsnr=2, bandtype=Btype, degamp=args.degamp, degphase=args.degphase,
                    caltable=bc, solint='inf', refant=args.refant, timerange=args.timerange)
            plotcal(bc, xaxis='chan', yaxis='amp', figfile="{}.amp.png".format(bc), showgui=False)
            plotcal(bc, xaxis='chan', yaxis='phase', figfile="{}.phs.png".format(bc), showgui=False)

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
                print("...Saving bandpass to {}.csv".format(bc))
            else:
                print("couldn't access bandpass data. Note BPOLY solutions not currently compatible w/ caltable2calfits.py")

            applycal(KGsplit, gaintable=[bc])
            Bsplit = os.path.join(out_dir, "{}.Bsplit".format(base_ms))
            files = glob.glob("{}*".format(Bsplit))
            for f in files:
                if os.path.exists(f):
                    try:
                        shutil.rmtree(f)
                    except OSError:
                        os.remove(f)

            split(KGsplit, Bsplit, datacolumn="corrected")
            fixvis(Bsplit, Bsplit, phasecenter=fixdir)
            ft(Bsplit, complist="{}.cl".format(args.source), usescratch=True)

    # clean
    im_stem = os.path.join(out_dir, base_ms + '.' + args.source + source_ext)
    print("...performing clean for {}.* output files".format(im_stem))
    source_files = glob.glob(im_stem+'*')
    if len(source_files) > 0:
        for sf in source_files:
            if os.path.exists(sf):
                try:
                    shutil.rmtree(sf)
                except OSError:
                    os.remove(sf)

    # create mfs image
    print("...running MFS clean")
    if args.cleanspw is None:
        cleanspw = "0:100~924"
    else:
        cleanspw = args.cleanspw
    clean(vis=Bsplit, imagename=im_stem, spw=cleanspw, niter=args.niter, weighting='briggs', robust=-0.5, imsize=[args.imsize, args.imsize],
          cell=['{}arcsec'.format(args.pxsize)], mode='mfs', timerange=args.timerange)
    exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))

    # create spectrum
    if args.spec_cube:
        print("...running spectral cube clean")
        dchan = args.spec_dchan
        for i, chan in enumerate(np.arange(100, 924, dchan)):
            clean(vis=Bsplit, imagename=im_stem+'.spec{}'.format(chan), niter=0, spw="0:{}~{}".format(chan, chan+dchan-1),
                    weighting='briggs', robust=0, imsize=[args.imsize, args.imsize], timerange=args.timerange,
                    cell=['{}arcsec'.format(args.pxsize)], mode='mfs')#, mask='circle[[{}h{}m{}s, {}d{}m{}s ], 7deg]'.format(*(ra+dec)))
            exportfits(imagename='{}.spec{}.image'.format(im_stem, chan), fitsimage='{}.spec{}.fits'.format(im_stem, chan))

    # make uvdist plot
    if args.plot_uvdist:
        print("...plotting uvdistance")
        # load visibility amplitudes
        ms.open(Bsplit)
        data = ms.getdata(["amplitude", "antenna1", "antenna2", "uvdist", "axis_info"], ifraxis=True)
        amps = data['amplitude']
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
        # average across time
        amps = np.nanmedian(amps, 2)
        uvdist = np.nanmedian(uvdist, 1)
        # average into channel bins
        Nchanbin = 6
        chan_bins = np.digitize(np.arange(1024), np.linspace(0, 1024, Nchanbin, endpoint=True))
        freq_bins = np.linspace(100., 200., Nchanbin, endpoint=True)
        freq_bins = (freq_bins + np.nanmedian(np.diff(freq_bins))/2)[:-1]
        amps_temp = []
        for i in range(1, Nchanbin):
            select = np.where(chan_bins == i)[0]
            amps_temp.append(np.nanmedian(amps[select, :], axis=0))
        amps = np.array(amps_temp)

        # average into uvdist bins
        Nuvdbins = 21
        uvdist_bins = np.linspace(uvdist.min(), uvdist.max(), Nuvdbins, endpoint=True)
        uvd_bins = np.digitize(uvdist, uvdist_bins)
        uvdist_bins = np.around((uvdist_bins + np.nanmedian(np.diff(uvdist_bins))/2)[:-1], 1)
        amps_temp = []
        uvdist_temp = []
        for i in range(1, Nuvdbins):
            select = np.where(uvd_bins == i)[0]
            amps_temp.append(np.nanmedian(amps[:, select], axis=1))
            uvdist_temp.append(np.nanmedian(uvdist[select]))
        amps = np.array(amps_temp).T
        uvdist = np.array(uvdist_temp)

        # normalize amplitudes
        amps /= np.nanmedian(amps, axis=1).reshape(-1, 1)
        amps += np.linspace(0, 1, amps.shape[0]).reshape(-1, 1)

        # plot
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(10,7))
        ax.grid(True)
        p = ax.plot(uvdist_bins, amps.T, marker='s', markersize=6, alpha=0.8)
        ax.set_xlabel("baseline length (meters)", fontsize=16)
        ax.set_ylabel("amplitude (arbitrary units)", fontsize=16)
        ax.set_title("uvdist for {}".format(os.path.basename(msin)))
        ax.legend(p, map(lambda x: str(round(x, 1))+' MHz', freq_bins))
        ax.set_ylim(0, 3.0)
        fig.savefig(os.path.join(out_dir, "{}.uvdist.png".format(os.path.basename(msin))), bbox_inches='tight', pad=0.05)
        plt.close()



