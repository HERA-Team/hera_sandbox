#!/usr/bin/env python2.7
"""
abscal_pipeline.py
-----------------

Pipeline for absolute calibration
of HERA data.

Nick Kern
Jan 2018
"""
# Import Modules
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import hera_cal as hc
import hera_qm as hq
import uvtools as uvt
from pyuvdata import UVCal, UVData
import os, sys
import fnmatch
import cPickle
import astropy.stats as astats
import pygsm
import datetime
import glob
from memory_profiler import memory_usage
import omnical
import copy
import time
import json
import subprocess
import pathos
import aipy
from collections import OrderedDict
from astropy.time import Time
from scipy import interpolate
from source2file import source2file
from sklearn import linear_model
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
import aplpy

# import parameter file
from abscal_params import *
if overwrite:
    overwrite = '--overwrite'
else:
    overwrite = ''

# get source info
source_ra, source_dec = np.loadtxt("{}.loc".format(source), dtype=str)
source_ra, source_dec= source_ra.split(':'), source_dec.split(':')
source_ra = (float(source_ra[0]) + float(source_ra[1])/60. + float(source_ra[2])/3600.) * 15
source_dec = (float(source_dec[0]) + np.sign(float(source_dec[0]))*float(source_dec[1])/60. + np.sign(float(source_dec[0]))*float(source_dec[2])/3600.)

def echo(message, type=0):
    if type == 0:
        print message
    elif type == 1:
         print "\n"+message+"\n"+"-"*60

# get files
uv_files = sorted(glob.glob("{}/zen*.HH.uv".format(data_path)))
xx_files = sorted([x for x in uv_files if 'xx' in x])
xx_bases = map(lambda x: os.path.basename(x), xx_files)
yy_files = sorted([x for x in uv_files if 'yy' in x])
yy_bases = map(lambda x: os.path.basename(x), yy_files)
jd_files = map(lambda x: '.'.join(os.path.basename(x).split('/')[-1].split('.')[1:3]), xx_files)
Nfiles = len(xx_files)
devnull = open(os.devnull, 'w')

# open output files
abs_out = open("{}/abs_out.txt".format(data_path), 'w')
abs_err = open("{}/abs_err.txt".format(data_path), 'w')

if run_abscal:
    def abscal(uv_files, pol=-5):
        pol2str = {-5:'xx', -6:'yy'}
        polstr = pol2str[pol]

        # echo info
        echo("Running abscal on field {} in JD {} for pol {}".format(field, JD, pol), type=1)

        # get source info
        echo("getting source info", type=1)
        (lst, jd, utc_range, utc_center, source_files, source_utc_range) = source2file(source_ra, duration=duration, offset=0.0, start_jd=JD,
                                                                                       jd_files=uv_files, get_filetimes=True, verbose=True)
        utc_center = "'" + '/'.join(utc_center.split('/')[:-1]) + ' ' + utc_center.split('/')[-1] + "'"
        source_files = list(source_files)
        Nsf = len(source_files)

        # make flux model
        echo("making flux model", type=1)
        cmd = "casa --nologger --nocrashreport --nogui --agg -c {} --image --freqs 100,200,1024 --cell 150arcsec --imsize 512".format(complist)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("flux model exit {}".format(out))

        # apply PB correction to flux model
        echo("applying PB to flux model")
        cmd = "pbcorr.py --beamfile {} --outdir ./ --pol {} --time {} --ext pbcorr --lon {} --lat {} {} --multiply {}.cl.fits"\
              "".format(beamfile, pol, utc_center, lon, lat, overwrite, field)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("flux model PB corr exit {}".format(out))

        # import to CASA image
        echo("importing model file to CASA image format")
        cmd = '''casa --nologger --nocrashreport --nogui --agg -c ''' \
              '''"importfits('{}.cl.pbcorr.fits', '{}.cl.pbcorr.image', overwrite=True)"'''.format(field, field)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("flux model CASA import exit {}".format(out))
        modelfile = field + '.cl.pbcorr.image'

        if rfi_flag:
            echo("RFI flagging data", type=1)
            for i, sf in enumerate(source_files):
                cmd = "xrfi_run.py --algorithm 'xrfi' --kf_size 13 --summary {}".format(sf)
                out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
                echo("RFI flag on {} exit {}".format(sf, out))

            echo("applying RFI flags")
            for i, sf in enumerate(source_files):
                cmd = "xrfi_apply.py --ext R {} --flag_file {}.flags.npz {}".format(overwrite, sf, sf)
                out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
                echo("RFI apply on {} exit {}".format(sf, out))

            # update source file(s) paths
            for i in range(Nsf): source_files[i] = source_files[i] + 'R'

        # convert to uvfits
        echo("converting miriad to uvfits", type=1)
        for i, sf in enumerate(source_files):
            cmd = "miriad_to_uvfits.py {} {}".format(overwrite, sf)
            out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
            echo("miriad_to_uvfits on {} exit {}".format(sf, out))            
        source_uvfits = map(lambda x: x+'.uvfits', source_files)

        # combine uvfits
        if len(source_uvfits) > 1:
            echo("combining uvfits", type=1)
            cmd = "combine_uvfits.py {} {}".format(overwrite, ' '.join(source_uvfits))
            out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
            echo("combine uvfits exit {}".format(out))
        source_uvfile = source_files[0]
        source_uvfits = source_uvfits[0]

        # absolute calibration
        echo("absolutely calibrating", type=1)
        cmd = "casa --nologger --nocrashreport --nogui --agg -c sky_image.py " \
              "--msin {} --source {} --model_im {} --refant {} {} --imsize {} --pxsize {} --niter {} --timerange {} " \
              "--KGcal --Acal --BPcal --image_mfs --image_model --plot_uvdist" \
              "".format(source_uvfits, source, modelfile, refant, ex_ants, imsize, pxsize, niter, source_utc_range)
        if casa_rflag:
            cmd += ' --rflag'
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("abscal exit {}".format(out))

        # convert to calfits
        echo("converting CASA caltables to calfits", type=1)
        calfits_fname = source_uvfile + '.ms.abs.calfits'
        dly_file = source_uvfile + '.ms.K.cal.npz'
        phs_file = source_uvfile + '.ms.Gphs.cal.npz'
        amp_file = source_uvfile + '.ms.Gamp.cal.npz'
        bp_file = source_uvfile + '.ms.B.cal.npz'
        cmd = "skynpz2calfits.py --fname {} --uv_file {} --dly_file {} --phs_file {} --amp_file {} --bp_file {} " \
              "--plot_bp --plot_phs --plot_amp --plot_dlys --bp_medfilt --medfilt_kernel 13 {} " \
              "--bp_gp_max_dly {} --bp_gp_thin 4 {}".format(calfits_fname, source_uvfile, dly_file, phs_file, amp_file, 
                                                            bp_file, bp_gp_smooth, bp_gp_mx_dly, overwrite)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("skynpz2calfits exit {}".format(out))

        # primary beam correct MFS image
        echo("pb correct mfs image", type=1)
        cmd = "pbcorr.py --beamfile {} --pol {} --time {} --ext pbcorr --lon {} --lat {} --spec_cube {} {}" \
              "".format(beamfile, pol, utc_center, lon, lat, overwrite, source_uvfile+".ms.{}.fits".format(source))
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("pb correction exit {}".format(out))

        # plot source and model source
        source_im = source_uvfile + ".ms.{}".format(source)
        model_im = source_uvfile + ".ms.model.{}".format(source)
        vmin = -5
        vmax = 40

        fig = plt.figure(figsize=(14,10))

        f1 = aplpy.FITSFigure(source_im+".fits", figure=fig, subplot=[0.0, 0.55, 0.45, 0.45])
        f1.show_colorscale(cmap='nipy_spectral', vmin=vmin, vmax=vmax)
        f1.add_grid()
        f1.add_beam()
        f1.recenter(source_ra, source_dec, width=15, height=15)
        f1.add_colorbar()
        f1.colorbar.set_axis_label_text("Jy/beam")
        f1.set_title("{} {} MFS Image".format(source, polstr))

        f2 = aplpy.FITSFigure(source_im+".pb.fits", figure=fig, subplot=[0.55, 0.55, 0.45, 0.45])
        f2.show_colorscale(cmap='nipy_spectral', vmin=0, vmax=1)
        f2.add_grid()
        f2.recenter(source_ra, source_dec, width=15, height=15)
        f2.add_colorbar()
        f2.colorbar.set_axis_label_text("PB Beam Response")
        f2.set_title("Primary Beam")

        f3 = aplpy.FITSFigure(source_im+".pbcorr.fits", figure=fig, subplot=[0.0, 0.0, 0.45, 0.45])
        f3.show_colorscale(cmap='nipy_spectral', vmin=vmin, vmax=vmax)
        f3.add_grid()
        f3.add_beam()
        f3.recenter(source_ra, source_dec, width=15, height=15)
        f3.add_colorbar()
        f3.colorbar.set_axis_label_text("Jy/beam")
        f3.set_title("{} {} MFS Image + PB Correction".format(source, polstr))

        f4 = aplpy.FITSFigure(model_im+".fits", figure=fig, subplot=[0.55, 0.0, 0.45, 0.45])
        f4.show_colorscale(cmap='nipy_spectral', vmin=vmin, vmax=vmax)
        f4.add_grid()
        f4.add_beam()
        f4.recenter(source_ra, source_dec, width=15, height=15)
        f4.add_colorbar()
        f4.colorbar.set_axis_label_text("Jy/beam")
        f4.set_title("Model MFS Image")

        fig.savefig("{}/{}_{}_MFS.png".format(data_path, source, polstr), dpi=150, bbox_inches='tight')
        plt.close()

        return calfits_fname

    abs_xx_calfits = abscal(xx_files, pol=-5)
    abs_yy_calfits = abscal(yy_files, pol=-6)

# apply abscal
if apply_abscal:
    def apply_abs(uv_file, calfits, pol='xx'):
        cmd = "omni_apply.py -p {} --omnipath {} --extension X {} {}".format(pol, calfits, overwrite, uv_file)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("omni apply on {} exit {}".format(uv_file, out))
        return out
    
    # assign mp pooling
    if multiprocess is True:
        pool = pathos.multiprocessing.Pool(Nproc)
        M = pool.map
    else:
        M = map

    echo("applying abscal calfits on XX pol", type=1)
    abs_apply_xx_exits = M(lambda x: apply_abs(x, abs_xx_calfits, 'xx'), xx_files)
    echo("applying abscal calfits on YY pol", type=1)
    abs_apply_yy_exits = M(lambda x: apply_abs(x, abs_yy_calfits, 'yy'), yy_files)

# source spectrum
if source_spectrum:
    # get flux calibrated files
    xx_files = map(lambda x: x + 'X', xx_files)
    yy_files = map(lambda x: x + 'X', yy_files)

    def source_extract(uv_files, pol=-5):
        pol2str = {-5:'xx', -6:'yy'}
        polstr = pol2str[pol]

        # get source info
        echo("getting source info", type=1)
        (lst, jd, utc_range, utc_center, source_files, source_utc_range) = source2file(source_ra, duration=duration, offset=0.0, start_jd=JD,
                                                                                       jd_files=uv_files, get_filetimes=True, verbose=True)
        utc_center = "'" + '/'.join(utc_center.split('/')[:-1]) + ' ' + utc_center.split('/')[-1] + "'"
        source_files = list(source_files)
        Nsf = len(source_files)

        # make flux model
        echo("making flux model", type=1)
        cmd = "casa --nologger --nocrashreport --nogui --agg -c {} --image --freqs 100,200,1024 --cell 150arcsec --imsize 512".format(complist)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("flux model exit {}".format(out))

        # import to CASA image
        echo("importing model file to CASA image format")
        cmd = '''casa --nologger --nocrashreport --nogui --agg -c ''' \
              '''"importfits('{}.cl.fits', '{}.cl.image', overwrite=True)"'''.format(field, field)
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("flux model CASA import exit {}".format(out))
        modelfile = field + '.cl.image'

        if spec_rfi_flag:
            echo("RFI flagging data", type=1)
            for i, sf in enumerate(source_files):
                cmd = "xrfi_run.py --algorithm 'xrfi' --kf_size 13 --summary {}".format(sf)
                out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
                echo("RFI flag on {} exit {}".format(sf, out))

            echo("applying RFI flags")
            for i, sf in enumerate(source_files):
                cmd = "xrfi_apply.py --ext R {} --flag_file {}.flags.npz {}".format(overwrite, sf, sf)
                out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
                echo("RFI apply on {} exit {}".format(sf, out))

            # update source file(s) paths
            for i in range(Nsf): source_files[i] = source_files[i] + 'R'

        # convert to uvfits
        echo("converting miriad to uvfits", type=1)
        for i, sf in enumerate(source_files):
            cmd = "miriad_to_uvfits.py {} {}".format(overwrite, sf)
            out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
            echo("miriad_to_uvfits on {} exit {}".format(sf, out))
        source_uvfits = map(lambda x: x+'.uvfits', source_files)

        # combine uvfits
        if len(source_uvfits) > 1:
            echo("combining uvfits", type=1)
            cmd = "combine_uvfits.py {} {}".format(overwrite, ' '.join(source_uvfits))
            out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
            echo("combine uvfits exit {}".format(out))
        source_uvfile = source_files[0]
        source_uvfits = source_uvfits[0]

        # spectral cube imaging
        echo("imaging spectral cube", type=1)
        cmd = "casa --nologger --nocrashreport --nogui --agg -c sky_image.py " \
              "--msin {} --source {} --model_im {} {} --imsize {} --pxsize {} --niter {} --timerange {} " \
              "--image_mfs --image_model --spec_cube --spec_dchan {}" \
              "".format(source_uvfits, source, modelfile, ex_ants, imsize, pxsize, niter, source_utc_range, spec_dchan)
        if spec_casa_rflag:
            cmd += ' --rflag'
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("spec_cube exit {}".format(out))
        source_im = source_uvfile + ".ms.{}".format(source)
        model_im = source_uvfile + ".ms.model.{}".format(source)

        # primary beam correction
        echo("applying primary beam correction to spectral cube", type=1)
        cmd = "pbcorr.py --beamfile {} --pol {} --time {} --ext pbcorr --lon {} --lat {} --spec_cube {} {}" \
              "".format(beamfile, pol, utc_center, lon, lat, overwrite, source_im + ".spec????.fits")
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("pb correction exit {}".format(out))

        # source extraction
        echo("running source extraction on data", type=1)
        cmd = "source_extract.py --source {} --radius 1 {} --gaussfit_mult {} --ext .spectrum.tab --rms_max_r {} --rms_min_r {} {}" \
              "".format(source, overwrite, gf_mult, rms_maxr, rms_minr, source_im + ".spec????.pbcorr.fits")
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("source extract exit {}".format(out))

        echo("running source extraction on model", type=1)
        cmd = "source_extract.py --source {} --radius 1 {} --gaussfit_mult {} --ext .spectrum.tab --rms_max_r {} --rms_min_r {} {}" \
                "".format(source, overwrite, gf_mult, rms_maxr, rms_minr, model_im + ".spec????.fits")
        out = subprocess.call(cmd, shell=True, stdout=abs_out, stderr=abs_err)
        echo("model extract exit {}".format(out))

        # plot spectrum
        echo("plotting spectrum", type=1)
        # load spectrum # freqs (MHz)  peak flux (Jy)  flux err    peak gaussfit (Jy)  integ gaussfit (Jy)
        source_spec = np.loadtxt("{}.spectrum.tab".format(source_im))
        source_freqs = source_spec[:, 0]
        source_flux = source_spec[:, 3]
        source_flux_err = source_spec[:, 2]

        # load source model spectrum
        source_model_spec = np.loadtxt("{}.spectrum.tab".format(model_im))
        source_model_freqs = source_model_spec[:, 0]
        source_model_flux = source_model_spec[:, 3]

        # fit power law to model spectrum
        model = make_pipeline(PolynomialFeatures(1), linear_model.RANSACRegressor())
        model.fit(np.log10(source_model_freqs).reshape(-1, 1), np.log10(source_model_flux))
        source_model_flux_fit = 10**(model.predict(np.log10(source_model_freqs).reshape(-1, 1)))

        # plot spectra
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax.grid()
        ax.set_xlabel("Frequency [MHz]", fontsize=18)
        ax.set_ylabel("Flux [Jansky]", fontsize=18)
        ax.set_title("{} Recovered {} Spectrum after Absolute Calibration".format(source, polstr), fontsize=14)
        _ = [tl.set_size(14) for tl in ax.get_xticklabels()]
        _ = [tl.set_size(14) for tl in ax.get_yticklabels()]

        p1, = ax.plot(source_model_freqs, source_model_flux, color='darkorange', ls='', marker='o', alpha=0.75)
        p2, = ax.plot(source_model_freqs, source_model_flux_fit, color='grey', lw=3)
        p0 = ax.errorbar(source_freqs, source_flux, yerr=source_flux_err, fmt='o', color='steelblue', ms=4)

        ax.legend([p0, p1, p2], ["Recovered Data Flux", "Recovered Model Flux", "Recovered Model Flux Fit"], fontsize=20)
        ax.set_ylim(0, np.max(np.concatenate([source_model_flux, source_flux, source_model_flux_fit])*1.8))

        fig.savefig("{}/{}_{}_spec.png".format(data_path, source, polstr), dpi=150, bbox_inches='tight')
        echo("saving {}_{}_spec.png".format(source, polstr))
        plt.close()

    source_extract(xx_files, pol=-5)
    source_extract(yy_files, pol=-6)


sys.exit(0)
