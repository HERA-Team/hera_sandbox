#!/usr/bin/env python2.7
"""
redcal_pipeline.py
-------------

Reprocessing for H1C.
Requires a redcal_params.py 
files in working directory.

Nicholas Kern
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
import argparse

# Import Params
from redcal_params import *

sys.exit(0)
# get info
uvd = UVData()
uvd.read_miriad(xx_files[0])
ants = uvd.antenna_numbers.tolist()
Nants = len(ants)
freqs = uvd.freq_array.squeeze()
Nfreqs = len(freqs)
HHaa = hc.utils.get_aa_from_uv(uvd)
info = hc.omni.aa_to_info(HHaa)
antpos, ants = uvd.get_ENU_antpos()
red_bls = info.get_reds()

####################
### Run FirstCal ###
####################

## Load Ant Mets xants ##
xants_xx = np.loadtxt(os.path.join(out_dir, 'xants_xx.csv'), dtype=np.int, delimiter=',')
xants_yy = np.loadtxt(os.path.join(out_dir, 'xants_yy.csv'), dtype=np.int, delimiter=',')

xants_xx = [','.join(xants_xx.astype(np.str)) for i in range(Nfiles)]
xants_yy = [','.join(xants_yy.astype(np.str)) for i in range(Nfiles)]

fc_xx_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.first.calfits'), xx_files)
fc_yy_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.first.calfits'), yy_files)

if run_firstcal is True:

    fc_out = open(os.path.join(out_dir, "fc_out.txt"), 'w')
    fc_err = open(os.path.join(out_dir, "fc_err.txt"), 'w')
    def fc_run(inp, Nfiles=Nfiles, out_dir=out_dir, fc_out=fc_out, fc_err=fc_err):
        # check if file exists
        i = inp[0]
        uvfile = inp[1][0]
        xants = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        # check if file exists
        uvbase = os.path.basename(uvfile)
        fc_fname = uvbase + ".first.calfits"
        if os.path.exists(os.path.join(out_dir, fc_fname)) and overwrite is False:
            return 0
        pol = os.path.basename(uvfile).split('.')[3]
        cmd = "firstcal_run.py -p {} --reds_tolerance=5.0 --maxiter=3 --ex_ants={} --overwrite --outpath={} {}".format(pol, xants, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=fc_out, stderr=fc_err)
        return out

    # run fc
    print_message("running firstcal on XX pol", type=1)
    fc_xx_exits = M(fc_run, enumerate(zip(xx_files, xants_xx)[:1]))
    print_message("running firstcal on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(zip(yy_files, xants_yy)[:1]))

    fc_out.close()
    fc_err.close()
    print "fc_xx_exits:", fc_xx_exits
    print "fc_yy_exits:", fc_yy_exits

fcm_xx_files = map(lambda x: x + '.firstcal_metrics.json', fc_xx_files)
fcm_yy_files = map(lambda x: x + '.firstcal_metrics.json', fc_yy_files)

if run_fcmets is True:
    print_message("running firstcal metrics", type=1)

    fc_out = open(os.path.join(out_dir, "fc_out.txt"), 'w')
    fc_err = open(os.path.join(out_dir, "fc_err.txt"), 'w')
    def fc_run(inp, Nfiles=Nfiles, out_dir=out_dir, fc_out=fc_out, fc_err=fc_err):
        i = inp[0]
        calfits_file = inp[1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "firstcal_metrics_run.py --metrics_path={} {}".format(out_dir, calfits_file)
        out = subprocess.call(cmd, shell=True, stdout=fc_out, stderr=fc_err)
        return out

    # run fc
    print_message("running firstcal metrics on XX pol", type=1)
    fc_xx_exits = M(fc_run, enumerate(fc_xx_files[:1]))
    print_message("running firstcal metrics on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(fc_yy_files[:1]))

    fc_out.close()
    fc_err.close()
    print "fc_xx_exits:", fc_xx_exits
    print "fc_yy_exits:", fc_yy_exits

if collate_fcmets is True:
    print_message("collating firstcal and firstcal metrics", type=1)

    # gather fc metrics
    fc_xx_xants = []
    fc_yy_xants = []
    fc_xx_rerun = []
    fc_yy_rerun = []

    mets_xx = OrderedDict({'file_times':[],
                            'ant_dly_std':OrderedDict(map(lambda a: (a, []), ants)),
                            'ant_dlys':OrderedDict(map(lambda a: (a, []), ants)),
                            'rot_ants':[]})

    mets_yy = OrderedDict({'file_times':[],
                            'ant_dly_std':OrderedDict(map(lambda a: (a, []), ants)),
                            'ant_dlys':OrderedDict(map(lambda a: (a, []), ants)),
                            'rot_ants':[]})

    for i, (fcxx, fcmxx) in enumerate(zip(fc_xx_files, fcm_xx_files)):
        xants = []
        if os.path.isfile(fcmxx) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fcxx)
            mets_xx['file_times'].extend(FC.times)
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmxx)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
                mets_xx['rot_ants'].append(FC.metrics['rot_ants'])
            fc_xx_xants.append(xants)
            if len(xants) > 0:
                fc_xx_rerun.append(i)
            for a in ants:
                if a in FC.ants:
                    mets_xx['ant_dly_std'][a].append(FC.metrics['ant_std'][a])
                    mets_xx['ant_dlys'][a].extend(FC.delays[FC.ants.tolist().index(a)])
                else:
                    mets_xx['ant_dly_std'][a].append(np.nan)
                    mets_xx['ant_dlys'][a].extend([np.nan for i in range(len(FC.times))])
        else:
            fc_xx_rerun.append(i)
            fc_xx_xants.append(xants)

    for i, (fcyy, fcmyy) in enumerate(zip(fc_yy_files, fcm_yy_files)):
        xants = []
        if os.path.isfile(fcmyy) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fcyy)
            mets_yy['file_times'].extend(FC.times)
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmyy)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
                mets_yy['rot_ants'].append(FC.metrics['rot_ants'])
            fc_yy_xants.append(xants)
            if len(xants) > 0:
                fc_yy_rerun.append(i)
            for a in ants:
                if a in FC.ants:
                    mets_yy['ant_dly_std'][a].append(FC.metrics['ant_std'][a])
                    mets_yy['ant_dlys'][a].extend(FC.delays[FC.ants.tolist().index(a)])
                else:
                    mets_yy['ant_dly_std'][a].append(np.nan)
                    mets_yy['ant_dlys'][a].extend([np.nan for i in range(len(FC.times))])
        else:
            fc_yy_rerun.append(i)
            fc_yy_xants.append(xants)


    # plot dly std
    def plot_dly_std(fname, mets):
        fig = plt.figure(figsize=(14,8), dpi=100)
        ax = fig.add_subplot(111)
        ax.grid(True)
        # plotting
        cm_func = plt.get_cmap('nipy_spectral')
        cm = cm_func(np.linspace(0, 0.95, len(ants)))
        p = map(lambda x: ax.plot(np.array(jd_files, np.float), x[1], marker='o', ls='-', c=cm[x[0]]), enumerate(np.array(mets['ant_dly_std'].values())))
        # axes
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        [t.set_rotation(20) for t in ax.get_xticklabels()]
        [t.set_fontsize(12) for t in ax.get_xticklabels()]
        [t.set_fontsize(12) for t in ax.get_yticklabels()]
        ax.set_xlabel('decimal of JD {}'.format(JD), fontsize=14)
        ax.set_ylabel('standard deviation [nanosec]', fontsize=16)
        ax.set_title("standard deviation of delay fluctuation per file", fontsize=14)
        ax1 = fig.add_axes([.95,0.1,0.05,0.8])
        ax1.axis('off')
        leg = ax1.legend(np.concatenate(p), ants, ncol=2)
        fig.savefig(fname)
        plt.close()

    plot_dly_std('dlystd_xx_{}.png'.format(JD), mets_xx)
    plot_dly_std('dlystd_yy_{}.png'.format(JD), mets_yy)

    # plot delay fluctuations
    def plot_dly_fluc(fname, mets):
        fig, ax = plt.subplots(1, 1, figsize=(14, 6), dpi=100)
        ax.grid(True)
        cm_func = plt.get_cmap('nipy_spectral')
        cm = cm_func(np.linspace(0, 0.95, len(ants)))
        p = map(lambda x: ax.plot(mets['file_times'], np.array(x[1])-np.nanmedian(x[1]), marker='.', ls='-', c=cm[x[0]]), enumerate(np.array(mets['ant_dlys'].values())))
        # axes
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        [t.set_rotation(20) for t in ax.get_xticklabels()]
        [t.set_fontsize(12) for t in ax.get_xticklabels()]
        [t.set_fontsize(12) for t in ax.get_yticklabels()]
        ax.set_xlabel('decimal of JD {}'.format(JD), fontsize=14)
        ax.set_ylabel('delay fluctuation [nanosec]', fontsize=16)
        ax.set_title("delay fluctuations across the observation", fontsize=14)
        ax1 = fig.add_axes([.95,0.1,0.05,0.8])
        ax1.axis('off')
        leg = ax1.legend(np.concatenate(p), ants, ncol=2)
        fig.savefig(fname)
        plt.close()

    plot_dly_fluc('dlyfluc_xx_{}.png'.format(JD), mets_xx)
    plot_dly_fluc('dlyfluc_yy_{}.png'.format(JD), mets_yy)


###################
### Run OmniCal ###
###################

# eliminate files
eliminate = False
if eliminate is True:
    select = np.where(np.array(jd_files, np.float) % JD  > 0.3)
    jd_files = np.array(jd_files)[select]
    xx_files = np.array(xx_files)[select]
    yy_files = np.array(yy_files)[select]
    xx_bases = np.array(xx_bases)[select]
    yy_bases = np.array(yy_bases)[select]
    fc_xx_files = np.array(fc_xx_files)[select]
    fc_yy_files = np.array(fc_yy_files)[select]
    xants_xx = np.array(xants_xx)[select]
    xants_yy = np.array(xants_yy)[select]

oc_xx_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.omni.calfits'), xx_files)
oc_yy_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.omni.calfits'), yy_files)

if run_omnical is True:
    print_message("running omnical", type=1)

    oc_out = open(os.path.join(out_dir, "oc_out.txt"), 'w')
    oc_err = open(os.path.join(out_dir, "oc_err.txt"), 'w')
    def oc_run(inp, Nfiles=Nfiles, out_dir=out_dir, oc_out=oc_out, oc_err=oc_err):
        i = inp[0]
        cffile = inp[1][0]
        uvfile = inp[1][1]
        xants = inp[1][2]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        # check if file exists
        uvbase = os.path.basename(uvfile)
        fc_fname = uvbase + ".omni.calfits"
        if os.path.exists(os.path.join(out_dir, fc_fname)) and overwrite is False:
            return 0
        pol = os.path.basename(uvfile).split('.')[3]
        cmd = "omni_run.py -p {} --ex_ants {} --overwrite --firstcal={} --omnipath={} {}".format(pol, xants, cffile, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=oc_out, stderr=oc_err)
        return out

    # run oc
    print_message("running omnical on XX pol", type=1)
    oc_xx_exits = M(oc_run, enumerate(zip(fc_xx_files, xx_files, xants_xx)[:1]))
    print_message("running omnical on YY pol", type=1)
    oc_yy_exits = M(oc_run, enumerate(zip(fc_yy_files, yy_files, xants_yy)[:1]))

    oc_out.close()
    oc_err.close()
    print "oc_xx_exits:", oc_xx_exits
    print "oc_yy_exits:", oc_yy_exits

ocm_xx_files = map(lambda x: x + '.omni_metrics.json', oc_xx_files)
ocm_yy_files = map(lambda x: x + '.omni_metrics.json', oc_yy_files)

if run_ocmets is True:
    print_message("running omnical metrics", type=1)

    oc_out = open(os.path.join(out_dir, "oc_out.txt"), 'w')
    oc_err = open(os.path.join(out_dir, "oc_err.txt"), 'w')
    def oc_run(inp, Nfiles=Nfiles, out_dir=out_dir, oc_out=oc_out, oc_err=oc_err):
        i = inp[0]
        fcfile = inp[1][0]
        ocfile = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "omnical_metrics_run.py --fc_files={} {}".format(fcfile, ocfile)
        out = subprocess.call(cmd, shell=True, stdout=oc_out, stderr=oc_err)
        return out

    # run ocm
    print_message("running omnical metrics on XX pol", type=1)
    oc_xx_exits = M(oc_run, enumerate(zip(fc_xx_files, oc_xx_files)[:1]))
    print_message("running omnical metrics on YY pol", type=1)
    oc_yy_exits = M(oc_run, enumerate(zip(fc_yy_files, oc_yy_files)[:1]))

    oc_out.close()
    oc_err.close()
    print "oc_xx_exits:", oc_xx_exits
    print "oc_yy_exits:", oc_yy_exits

oa_xx_files = map(lambda x: os.path.join(out_dir, x + 'O'), xx_bases)
oa_yy_files = map(lambda x: os.path.join(out_dir, x + 'O'), yy_bases)

if apply_omni is True:
    print_message("applying omni", type=1)

    oa_out = open(os.path.join(out_dir, "oa_out.txt"), 'w')
    oa_err = open(os.path.join(out_dir, "oa_err.txt"), 'w')
    def oa_run(inp, Nfiles=Nfiles, out_dir=out_dir, oa_out=oa_out, oa_err=oa_err):
        i = inp[0]
        ocfile = inp[1][0]
        uvfile = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "omni_apply.py -p=xx --omnipath={} --outpath={} --overwrite {} ".format(ocfile, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=oa_out, stderr=oa_err)
        return out

    # run omni apply
    print_message("running omnical apply on XX pol", type=1)
    oa_xx_exits = M(oa_run, enumerate(zip(oc_xx_files, xx_files)[:1]))
    print_message("running omnical apply on YY pol", type=1)
    oa_yy_exits = M(oa_run, enumerate(zip(oc_yy_files, yy_files)[:1]))

    oa_out.close()
    oa_err.close()
    print "oa_xx_exits:", oa_xx_exits
    print "oa_yy_exits:", oa_yy_exits

flag_xx_files = map(lambda x: os.path.join(out_dir, x+'.flags.npz'), xx_bases)
flag_yy_files = map(lambda x: os.path.join(out_dir, x+'.flags.npz'), yy_bases)

if run_rfi is True:
    print_message("running xrfi", type=1)

    rfi_out = open(os.path.join(out_dir, "rfi_out.txt"), 'w')
    rfi_err = open(os.path.join(out_dir, "rfi_err.txt"), 'w')
    def rfi_run(inp, Nfiles=Nfiles, out_dir=out_dir, rfi_out=rfi_out, rfi_err=rfi_err):
        i = inp[0]
        uvfile = inp[1][0]
        ex_ants = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "xrfi_run.py --ex_ants={} --xrfi_path={} {}".format(ex_ants, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=rfi_out, stderr=rfi_err)
        return out

    # run ocm
    print_message("running xrfi  on XX pol", type=1)
    rfi_xx_exits = M(rfi_run, enumerate(zip(xx_files, xants_xx)[:1]))
    print_message("running xrfi on YY pol", type=1)
    rfi_yy_exits = M(rfi_run, enumerate(zip(yy_files, xants_yy)[:1]))

    rfi_out.close()
    rfi_err.close()
    print "rfi_xx_exits:", rfi_xx_exits
    print "rfi_yy_exits:", rfi_yy_exits

or_xx_files = map(lambda x: x+"R", oa_xx_files)
or_yy_files = map(lambda x: x+"R", oa_yy_files)

if apply_rfi is True:
    print_message("applying RFI", type=1)

    rfi_out = open(os.path.join(out_dir, "rfi_out.txt"), 'w')
    rfi_err = open(os.path.join(out_dir, "rfi_err.txt"), 'w')
    def rfi_run(inp, Nfiles=Nfiles, out_dir=out_dir, rfi_out=rfi_out, rfi_err=rfi_err):
        i = inp[0]
        flagfile = inp[1][0]
        uvfile = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "xrfi_apply.py --xrfi_path={} --flag_file={} --overwrite {}".format(out_dir, flagfile, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=rfi_out, stderr=rfi_err)
        return out

    # run ocm
    print_message("running rfi apply on XX pol", type=1)
    rfi_xx_exits = M(rfi_run, enumerate(zip(flag_xx_files, oa_xx_files)[:1]))
    print_message("running rfi apply on YY pol", type=1)
    rfi_yy_exits = M(rfi_run, enumerate(zip(flag_yy_files, oa_yy_files)[:1]))

    rfi_out.close()
    rfi_err.close()
    print "rfi_xx_exits:", rfi_xx_exits
    print "rfi_yy_exits:", rfi_yy_exits


if plot_reds is True:
    print_message("plotting redundant visbilities", type=1)

    # get data
    def plot_vis(uvfile, pol, mode='amp', ants=ants, red_bls=red_bls):
        uvbase = os.path.basename(uvfile)
        file_jd = float('.'.join(uvbase.split('.')[1:3]))

        uvd  = UVData()
        uvd.read_miriad(uvfile)

        # get badants
        if pol is 'xx':
            badants = np.unique((','.join(xants_xx)).split(',')).astype(int)
        elif pol is 'yy':
            badants = np.unique((','.join(xants_yy)).split(',')).astype(int)

        # get baseline seps
        ants = list(ants)
        bs_seps = np.array(map(lambda x: antpos[ants.index(x[0][1])] - antpos[ants.index(x[0][0])], red_bls))
        
        # get smallest ew
        ew = np.where(np.abs(bs_seps[:, 1]) < 1)[0]
        small_ew = ew[np.where(np.abs(bs_seps[ew])[:, 0] == np.abs(bs_seps[ew])[:, 0].min())[0]][0]
        bl_group = red_bls[small_ew]
        Nbls = len(bl_group)

        Nside = 3
        Yside = int(np.ceil(float(Nbls+1)/Nside))

        fig, axes = plt.subplots(Yside, Nside, figsize=(14, 14*float(Yside)/Nside), dpi=75)
        fig.subplots_adjust(wspace=0.1, hspace=0.3)
        if mode == 'amp':
            fig.suptitle("Omnical Amplitude Waterfalls for {} Pol & JD = {}".format(pol, file_jd), fontsize=14)
        elif mode == 'phs':
            fig.suptitle("Omnical Phase Waterfalls for {} Pol & JD = {}".format(pol, file_jd), fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.95))

        vis_data = []
        k = 0
        for i in range(Yside):
            for j in range(Nside):
                ax = axes[i, j]
                if k < Nbls:
                    bl = bl_group[k]
                    data = uvd.get_data(bl)
                    if bl[0] in badants or bl[1] in badants:
                        pass
                    else:
                        vis_data.append(data)
                    if mode == 'amp':
                        ax.matshow(np.log10(np.abs(data)), vmin=-3, vmax=2, aspect='auto')
                    elif mode == 'phs':
                        ax.matshow(np.angle(data), vmin=-np.pi, vmax=np.pi, aspect='auto')
                    rfi_flags = np.array(uvd.get_flags(bl).copy(), np.bool)
                    rfi_flags = np.ma.masked_where(~rfi_flags, rfi_flags)
                    ax.matshow(rfi_flags, cmap='bone_r', aspect='auto')
                    ax.xaxis.set_ticks_position('bottom')
                    ax.set_title("{}".format(bl_group[k]), fontsize=12, y=1.01)
                elif k == Nbls:
                    vis_data = np.array(vis_data)
                    mean_vis = np.median(vis_data, axis=0)
                    if mode == 'amp':
                        ax.matshow(np.log10(np.abs(mean_vis)), vmin=-3, vmax=2, aspect='auto')
                    elif mode == 'phs':
                        ax.matshow(np.angle(mean_vis), vmin=-np.pi, vmax=np.pi, aspect='auto')
                    ax.xaxis.set_ticks_position('bottom')
                    ax.set_title('baseline average', fontsize=12, y=1.01)
                else:
                    ax.axis('off')
                    
                if j != 0:
                    ax.set_yticklabels([])
                else:
                    [t.set_fontsize(10) for t in ax.get_yticklabels()]
                    ax.set_ylabel('time integrations', fontsize=10)
                if i != Yside-1:
                    ax.set_xticklabels([])
                else:
                    [t.set_fontsize(10) for t in ax.get_xticklabels()]
                    ax.set_xlabel('freq channel', fontsize=10)
                    
                k += 1 

        if mode == 'amp':
            fig.savefig("{}.EW_{}_AMPwaterfall.png".format(uvbase, pol))
        elif mode == 'phs':
            fig.savefig("{}.EW_{}_PHSwaterfall.png".format(uvbase, pol))
        plt.close()

    # iterate through files and pols
    for uvxx, uvyy in zip(or_xx_files, or_yy_files):
        plot_vis(uvxx, 'xx', mode='amp')
        plot_vis(uvxx, 'xx', mode='phs')
        plot_vis(uvyy, 'yy', mode='amp')
        plot_vis(uvyy, 'yy', mode='phs')

avg_vis = False
if avg_vis is True:
    print_message("averaging visibilities", type=1)

    uvfile = or_xx_files[0]

    uvd = UVData()
    uvd.read_miriad(uvfile)

    # create uvdata object
    uvd = UVData()

    uvd.Nants_telescope    = self.Nants
    uvd.Nants_data         = len(np.unique(bl_array))
    uvd.Nbls               = Nbls
    uvd.Ntimes             = Ntimes
    uvd.Nblts              = Nbls * Ntimes
    uvd.Nfreqs             = self.Nfreqs
    uvd.Npols              = self.Nxpols
    uvd.Nspws              = 1
    uvd.antenna_numbers    = self.ant_nums
    uvd.antenna_names      = np.array(self.ant_nums, dtype=str)
    uvd.ant_1_array        = np.concatenate([map(lambda x: x[0], bl_array) for i in range(Ntimes)])
    uvd.ant_2_array        = np.concatenate([map(lambda x: x[1], bl_array) for i in range(Ntimes)])
    uvd.baseline_array     = np.concatenate([np.arange(Nbls) for i in range(Ntimes)])
    uvd.freq_array         = (self.freqs * 1e6).reshape(1, -1)
    uvd.time_array         = np.repeat(np.array(JD_array)[:, np.newaxis], Nbls, axis=1).ravel()
    uvd.channel_width      = (self.freqs * 1e6)[1] - (self.freqs * 1e6)[0]
    uvd.data_array         = self.vis_data.reshape(uvd.Nblts, 1, self.Nfreqs, self.Nxpols)
    uvd.flag_array         = np.ones_like(uvd.data_array, dtype=np.bool)
    uvd.history            = " "
    uvd.instrument         = " "
    uvd.integration_time   = 10.7
    uvd.lst_array          = np.repeat(np.array(map(lambda x: self.JD2LST(self.loc, x), JD_array))[:, np.newaxis], Nbls, axis=1).ravel()
    uvd.nsample_array      = np.ones_like(uvd.data_array, dtype=np.float)
    uvd.object_name        = "GSM"
    uvd.phase_type         = "drift"
    uvd.polarization_array = np.array(map(lambda x: {"XX":-5,"YY":-6,"XY":-7,"YX":-8}[x], self.xpols))
    uvd.spw_array          = np.array([0])
    uvd.telescope_location = np.array([ 6378137.*np.cos(self.loc.lon)*np.cos(self.loc.lat),
                                        6378137.*np.sin(self.loc.lon)*np.cos(self.loc.lat),
                                        6378137.*np.sin(self.loc.lat) ])
    uvd.telescope_name     = " "
    uvd.uvw_array          = np.ones((uvd.Nblts, 3))
    uvd.vis_units          = "Jy"
    zen_dec, zen_ra        = self.loc.radec_of(0, np.pi/2)
    uvd.zenith_dec         = np.ones(uvd.Nblts) * zen_dec
    uvd.zenith_ra          = np.ones(uvd.Nblts) * zen_ra

    uvd.write_miriad(fname, clobber=clobber)





