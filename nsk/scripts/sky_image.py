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
a.add_argument('--msin', default=None, type=str, help='path to CASA measurement set. if fed a .uvfits, will convert to ms', required=True)
a.add_argument('--source', default=None, type=str, help='source name', required=True)
a.add_argument('--refant', default=None, type=str, help='reference antenna')
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument('--nocal', default=False, action='store_true', help='skip calibration and just make an image')
a.add_argument('--onlyKcal', default=False, action='store_true', help='only perform K calibration, then image')
a.add_argument('--onlyKGcal', default=False, action='store_true', help='only perform K & G calibration, then image')
a.add_argument('--source_ext', default=None, type=str, help="extension to source name in output image files")
a.add_argument('--spec_cube', default=False, action='store_true', help="image spectral cube as well as MFS.")
a.add_argument('--niter', default=50, type=int, help='number of clean iterations.')
a.add_argument('--pxsize', default=200, type=int, help='pixel (cell) scale in arcseconds')
a.add_argument('--imsize', default=350, type=int, help='number of pixels along a side of the output square image.')
a.add_argument('--bpoly', default=False, action='store_true', help="use BPOLY mode in bandpass")
a.add_argument('--degamp', default=4, type=int, help="amplitude polynomial degree for BPOLY")
a.add_argument('--degphase', default=1, type=int, help="phase polynomial degree for BPOLY")
a.add_argument('--noBPphase', default=False, action='store_true', help="set Bandpass phase solutions to identically zero")
args = a.parse_args()

if __name__ == "__main__":
    # get vars
    source = args.source
    refant = args.refant
    ex_ants = args.ex_ants
    rflag = args.rflag
    nocal = args.nocal
    onlyKcal = args.onlyKcal
    onlyKGcal = args.onlyKGcal
    spec_cube = args.spec_cube
    niter = args.niter
    source_ext = args.source_ext
    pxsize = args.pxsize
    imsize = args.imsize
    if source_ext is None:
        source_ext = ''

    # get paths
    msin = args.msin
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

    # fix vis
    print("...fix vis")
    ra, dec = np.loadtxt('{}.loc'.format(source), dtype=str)
    ra, dec = ra.split(':'), dec.split(':')
    fixvis(msin, msin, phasecenter='J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec)))

    # insert source model
    # set up calibrator source
    try:
        cl.addcomponent(flux=1, fluxunit='Jy', shape='point', dir='J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec)))
        cl.rename(os.path.join(out_dir, '{}.cl'.format(source)))
        cl.close()
    except:
        pass

    ft(msin, complist=os.path.join(out_dir, "{}.cl".format(source)), usescratch=True)

    # unflag
    if args.unflag is True:
        print("...unflagging")
        flagdata(msin, mode='unflag')

    # flag bad ants
    if ex_ants is not None:
        print("...flagging bad ants: {}".format(ex_ants))
        flagdata(msin, mode='manual', antenna=ex_ants)

    # flag autocorrs
    print("...flagging autocorrs")
    flagdata(msin, autocorr=True)

    # rflag
    if rflag is True:
        print("...rfi flagging")
        flagdata(msin, mode='rflag')

    # setup calibration tables
    kc = os.path.join(out_dir, base_ms+'.{}.cal'.format('K'))
    gc = os.path.join(out_dir, base_ms+'.{}.cal'.format('G'))

    if nocal is True:
        Bsplit = msin
    else:
        # perform initial K calibration (per-antenna delay)
        print("...performing K gaincal")
        if os.path.exists(kc):
            shutil.rmtree(kc)
        if os.path.exists("{}.png".format(kc)):
            os.remove("{}.png".format(kc))
        gaincal(msin, caltable=kc, gaintype="K", solint='inf', refant=refant, minsnr=3.0, spw="0:150~900", gaintable=[])
        plotcal(kc, xaxis='antenna', yaxis='delay', figfile='{}.png'.format(kc), showgui=False)

        # write delays to text file
        tb.open(kc)
        delays = tb.getcol('FPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]]
        delay_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
        tb.close()
        np.savetxt("{}.csv".format(kc), np.vstack([delay_ants, delays]).T, fmt="%12.8f", header="ant_num, delay [ns]", delimiter=", ")
        print("...Saving delays to {}.csv".format(kc))

        if onlyKcal:
            applycal(msin, gaintable=[kc])

        else:
            # perform initial G calibration (per-spw and per-pol gain)
            print("...performing G gaincal")
            if os.path.exists(gc):
                shutil.rmtree(gc)
            gaincal(msin, caltable=gc, gaintype='G', solint='inf', refant=refant, minsnr=3, calmode='p', spw="0:300~700", gaintable=[kc])
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
            applycal(msin, gaintable=[gc, kc])

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

        if onlyKGcal is True or onlyKcal is True:
            Bsplit = KGsplit
        else:
            # calibrated bandpass
            print("...performing B bandpass cal")
            bc = KGsplit+'.{}.cal'.format('B')
            Btype = "B"
            if args.bpoly:
                Btype="BPOLY"
            bandpass(vis=KGsplit, spw="0:100~900", minsnr=2, bandtype=Btype, degamp=args.degamp, degphase=args.degphase,
                    caltable=bc, solint='inf', refant=refant)
            plotcal(bc, xaxis='chan', yaxis='amp', figfile="{}.amp.png".format(bc), showgui=False)
            plotcal(bc, xaxis='chan', yaxis='phase', figfile="{}.phs.png".format(bc), showgui=False)

            # write bp to file
            try:
                tb.open(bc)
                flags = ~tb.getcol('FLAG')[0]
                flagged_ants = np.sum(flags, axis=0).astype(np.bool)
                flags = flags[:, flagged_ants]
                bp = tb.getcol('CPARAM')[0][:, flagged_ants]
                bp_ants = tb.getcol("ANTENNA1")[flagged_ants]
                tb.close()
                tb.open(bc+"/SPECTRAL_WINDOW")
                freqs = tb.getcol("CHAN_FREQ")
                tb.close()
                bp_real_data = bp.real
                bp_imag_data = bp.imag
                bp_data = np.concatenate([freqs/1e6, bp_real_data, bp_imag_data], axis=1)
                bp_ants_str = ', '.join(map(lambda x: str(x)+'r', bp_ants)) + ', ' + ', '.join(map(lambda x: str(x)+'i', bp_ants))
                np.savetxt("{}.csv".format(bc), bp_data, fmt="%12.8f", header="freq (MHz), {}".format(bp_ants_str), delimiter=", ")
                print("...Saving bandpass to {}.csv".format(gc))

            except:
                pass

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

    # clean
    im_stem = os.path.join(out_dir, base_ms + '.' + source + source_ext)
    print("...performing clean for {}.* output files".format(im_stem))
    source_files = glob.glob(im_stem+'*')
    if len(source_files) > 0:
        for sf in source_files:
            if os.path.exists(sf):
                try:
                    shutil.rmtree(sf)
                except OSError:
                    os.remove(sf)

    # craete mfs image
    print("...running MFS clean")
    clean(vis=Bsplit, imagename=im_stem, spw="0:200~850", niter=niter, weighting='briggs', robust=-0.5, imsize=[imsize, imsize],
          cell=['{}arcsec'.format(pxsize)], mode='mfs')
    exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))

    # create spectrum
    if spec_cube:
        print("...running spectral cube clean")
        dchan = 40
        for i, chan in enumerate(np.arange(100, 900, dchan)):
            clean(vis=Bsplit, imagename=im_stem+'.spec{}'.format(chan), niter=0, spw="0:{}~{}".format(chan, chan+dchan-1),
                    weighting='briggs', robust=0, imsize=[imsize, imsize],
                    cell=['{}arcsec'.format(pxsize)], mode='mfs')#, mask='circle[[{}h{}m{}s, {}d{}m{}s ], 7deg]'.format(*(ra+dec)))
            exportfits(imagename='{}.spec{}.image'.format(im_stem, chan), fitsimage='{}.spec{}.fits'.format(im_stem, chan))



