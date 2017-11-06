"""
reduction.py
-------------

Processing HERA Data Pipeline

Nicholas Kern
October, 2017
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

## set flags ##
T = True
F = False

run_antmets         = F
rerun_antmets       = F
collate_antmets     = F

run_firstcal        = F
run_fcmets          = F
add_fcxants         = F
rerun_firstcal      = F
rerun_fcmets        = F
collate_fcmets      = F

run_omnical         = F
run_ocmets          = F
apply_omni          = F

run_rfi             = F
apply_rfi           = F
plot_reds           = F

abs_cal             = F

multiprocess        = True
Nproc               = 10

# other vars
calfile = "hsa7458_v001"

def print_message(message, type=0):
    if type == 0:
        print message
    elif type == 1:
        print "\n"+message+"\n"+"-"*40

# get JD
JD = 2457680
try:
    JD = int(sys.argv[1])
except:
    pass

# assign mp pooling
if multiprocess is True:
    pool = pathos.multiprocessing.Pool(Nproc)
    M = pool.map
else:
    M = map

# Get data files
data_path = os.path.join("/data4/paper/HERA2015", str(JD))
out_dir = os.path.basename(data_path)
uv_files = sorted(glob.glob("{}/zen*.HH.uvc".format(data_path)))

uv_files = uv_files
xx_files = sorted([x for x in uv_files if 'xx' in x])
xx_bases = map(lambda x: os.path.basename(x), xx_files)
xy_files = sorted([x for x in uv_files if 'xy' in x])
yx_files = sorted([x for x in uv_files if 'yx' in x])
yy_files = sorted([x for x in uv_files if 'yy' in x])
yy_bases = map(lambda x: os.path.basename(x), yy_files)
jd_files = map(lambda x: '.'.join(os.path.basename(x).split('/')[-1].split('.')[1:3]), xx_files)

Nfiles = len(xx_files)
devnull = open(os.devnull, 'w')

# get info
uvd = UVData()
uvd.read_miriad(xx_files[0])
ants = np.unique(np.concatenate([uvd.ant_1_array, uvd.ant_2_array]))
Nants = len(ants)
freqs = uvd.freq_array.squeeze()
Nfreqs = len(freqs)
HHaa = hc.utils.get_aa_from_calfile(freqs, calfile)
info = hc.omni.aa_to_info(HHaa)
antpos = info.antloc
antpos = antpos[map(lambda x: info.ant_index(x), ants)]
red_bls = info.get_reds()

#######################
### Run Ant Metrics ###
#######################
am_out = open(os.path.join(out_dir, 'am_out.txt'), 'w')
am_err = open(os.path.join(out_dir, 'am_err.txt'), 'w')

def antmet_run(inp, Nfiles=Nfiles, out_dir=out_dir, am_out=am_out, am_err=am_err, calfile=calfile):
    i = inp[0]
    uvfile = inp[1]
    if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
    cmd = "ant_metrics_run.py -C {} -p xx,xy,yx,yy --metrics_path={} {}".format(calfile, out_dir, uvfile)
    out = subprocess.call(cmd, shell=True, stdout=am_out, stderr=am_err)
    return out

am_files = []
for xxf in xx_files:
    amf = os.path.basename(xxf).split('.')
    amf.pop(3)
    am_files.append(os.path.join(out_dir, '.'.join(amf))+'.ant_metrics.json')

if run_antmets is True:
    print_message("running ant metrics", type=1)
    # Run Ant Mets
    am_exits = M(antmet_run, enumerate(xx_files))
    print "am_exits:", am_exits

if rerun_antmets is True:
    rerun = np.arange(Nfiles)[map(lambda f: os.path.isfile(f) is False, am_files)]
    print_message("re-running ant metrics: {}".format(rerun), type=1)
    # Run Ant Mets
    am_exits = M(antmet_run, np.array(list(enumerate(xx_files)), np.object)[rerun])
    print "am_exits:", am_exits

am_out.close()
am_err.close()

if collate_antmets is True:
    # collect badants
    badxx = set()
    badyy = set()
    xants_xx = []
    xants_yy = []
    # loop over every file
    for i, amf in enumerate(am_files):
        if os.path.isfile(amf) is True:
            am = hq.ant_metrics.load_antenna_metrics(amf)
            xx = []
            yy = []
            for xant in am['xants']:
                if xant[1] == 'x':
                    xx.append(xant[0])
                elif xant[1] == 'y':
                    yy.append(xant[0])
            xants_xx.append([a in xx for a in ants])
            xants_yy.append([a in yy for a in ants])
            badxx.update(xx)
            badyy.update(yy)
        else:
            xants_xx.append(())
            xants_yy.append(())
    # for antmets that failed, assign all badants
    for i, o in enumerate(xants_xx):
        if o is ():
            xants_xx[i] = [a in ants for a in sorted(badxx)]
    for i, o in enumerate(xants_yy):
        if o is ():
            xants_yy[i] = [a in ants for a in sorted(badyy)]

    xants_xx = np.array(xants_xx)
    xants_yy = np.array(xants_yy)

    # write
    np.savetxt(os.path.join(out_dir, 'xants_xx.tab'), xants_xx, fmt="%d", delimiter='\t', header="ant_num")
    np.savetxt(os.path.join(out_dir, 'xants_yy.tab'), xants_yy, fmt="%d", delimiter='\t', header="ant_num")


####################
### Run FirstCal ###
####################

## Load Ant Mets xants ##
xants_xx = np.loadtxt(os.path.join(out_dir, 'xants_xx.tab'), dtype=np.bool)
xants_yy = np.loadtxt(os.path.join(out_dir, 'xants_yy.tab'), dtype=np.bool)

## Load Hard Coded xants ##
try:
    hard_xants_xx = np.loadtxt(os.path.join(out_dir, 'hard_xants_xx.tab'), dtype=np.int)
    hard_xants_yy = np.loadtxt(os.path.join(out_dir, 'hard_xants_yy.tab'), dtype=np.int)

    xants_xx += np.array([True if a in hard_xants_xx else False for a in ants])
    xants_yy += np.array([True if a in hard_xants_yy else False for a in ants])

except:
    pass

xants_xx = [','.join((ants[x]).astype(np.str)) for x in xants_xx]
xants_yy = [','.join((ants[x]).astype(np.str)) for x in xants_yy]

fc_xx_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.first.calfits'), xx_files)
fc_yy_files = map(lambda x: os.path.join(out_dir, os.path.basename(x)+'.first.calfits'), yy_files)

if run_firstcal is True:

    fc_out = open(os.path.join(out_dir, "fc_out.txt"), 'w')
    fc_err = open(os.path.join(out_dir, "fc_err.txt"), 'w')
    def fc_run(inp, Nfiles=Nfiles, out_dir=out_dir, fc_out=fc_out, fc_err=fc_err):
        i = inp[0]
        uvfile = inp[1][0]
        xants = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        pol = os.path.basename(uvfile).split('.')[3]
        cmd = "firstcal_run.py -C {} -p {} --ex_ants {} --overwrite --outpath={} {}".format(calfile, pol, xants, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=fc_out, stderr=fc_err)
        return out

    # run fc
    print_message("running firstcal on XX pol", type=1)
    fc_xx_exits = M(fc_run, enumerate(zip(xx_files, xants_xx)))
    print_message("running firstcal on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(zip(yy_files, xants_yy)))

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
    fc_xx_exits = M(fc_run, enumerate(fc_xx_files))
    print_message("running firstcal metrics on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(fc_yy_files))

    fc_out.close()
    fc_err.close()
    print "fc_xx_exits:", fc_xx_exits
    print "fc_yy_exits:", fc_yy_exits

    # gather fc metrics
    xants_xx = np.loadtxt(os.path.join(out_dir, 'xants_xx.tab'), dtype=np.bool)
    xants_yy = np.loadtxt(os.path.join(out_dir, 'xants_yy.tab'), dtype=np.bool)
    fc_xx_xants = []
    fc_yy_xants = []
    fc_xx_rerun = []
    fc_yy_rerun = []

    for i, fcmxx in enumerate(fcm_xx_files):
        xants = []
        if os.path.isfile(fcmxx) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fc_xx_files[i])
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmxx)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
            fc_xx_xants.append(xants)
            if len(xants) > 0:
                fc_xx_rerun.append(i)
        else:
            fc_xx_rerun.append(i)
            fc_xx_xants.append(xants)

    for i, fcmyy in enumerate(fcm_yy_files):
        xants = []
        if os.path.isfile(fcmyy) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fc_yy_files[i])
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmyy)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
            if len(xants) > 0:
                fc_yy_rerun.append(i)
            fc_yy_xants.append(xants)
        else:
            fc_yy_rerun.append(i)
            fc_yy_xants.append(xants)

    fc_xx_xants = np.array(map(lambda x: [a in x for a in ants], fc_xx_xants)) + xants_xx
    fc_yy_xants = np.array(map(lambda x: [a in x for a in ants], fc_yy_xants)) + xants_yy

    # write new xants to file
    if add_fcxants is True:
        np.savetxt(os.path.join(out_dir, 'xants_xx.tab'), fc_xx_xants, fmt="%d", delimiter='\t', header="ant_num")
        np.savetxt(os.path.join(out_dir, 'xants_yy.tab'), fc_yy_xants, fmt="%d", delimiter='\t', header="ant_num")

    # write reruns to file
    fc_xx_rerun_str = []
    for i in fc_xx_rerun:
        fc_xx_rerun_str.append( str(i) + '|' + xx_files[i] + '|' + ','.join(np.array(ants[fc_xx_xants[i]], np.str)) )
    fc_yy_rerun_str = []
    for i in fc_yy_rerun:
        fc_yy_rerun_str.append( str(i) + '|' + yy_files[i] + '|' + ','.join(np.array(ants[fc_yy_xants[i]], np.str)) )
    
    np.savetxt(os.path.join(out_dir, 'fc_xx_rerun.tab'), np.array(fc_xx_rerun_str), fmt="%s", header='index | uv file | xants')
    np.savetxt(os.path.join(out_dir, 'fc_yy_rerun.tab'), np.array(fc_yy_rerun_str), fmt="%s", header='index | uv file | xants')

xants_xx = np.loadtxt(os.path.join(out_dir, 'xants_xx.tab'), dtype=np.bool)
xants_yy = np.loadtxt(os.path.join(out_dir, 'xants_yy.tab'), dtype=np.bool)

## Load Hard Coded xants ##
try:
    hard_xants_xx = np.loadtxt(os.path.join(out_dir, 'hard_xants_xx.tab'), dtype=np.int)
    hard_xants_yy = np.loadtxt(os.path.join(out_dir, 'hard_xants_yy.tab'), dtype=np.int)

    xants_xx += np.array([True if a in hard_xants_xx else False for a in ants])
    xants_yy += np.array([True if a in hard_xants_yy else False for a in ants])

except:
    pass

xants_xx = [','.join((ants[x]).astype(np.str)) for x in xants_xx]
xants_yy = [','.join((ants[x]).astype(np.str)) for x in xants_yy]

if rerun_firstcal is True:
    print_message("re-running firstcal", type=1)

    fc_out = open(os.path.join(out_dir, "fc_out.txt"), 'w')
    fc_err = open(os.path.join(out_dir, "fc_err.txt"), 'w')
    def fc_run(inp, Nfiles=Nfiles, out_dir=out_dir, fc_out=fc_out, fc_err=fc_err):
        i = inp[0]
        uvfile = inp[1][0]
        xants = inp[1][1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        pol = os.path.basename(uvfile).split('.')[3]
        cmd = "firstcal_run.py -C {} -p {} --ex_ants {} --overwrite --outpath={} {}".format(calfile, pol, xants, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=fc_out, stderr=fc_err)
        return out

    # run fc
    fc_xx_rerun_str = np.loadtxt(os.path.join(out_dir, 'fc_xx_rerun.tab'), dtype=np.str)
    fc_xx_rerun_str = map(lambda x: x.split('|')[1:], fc_xx_rerun_str)
    fc_yy_rerun_str = np.loadtxt(os.path.join(out_dir, 'fc_yy_rerun.tab'), dtype=np.str)
    fc_yy_rerun_str = map(lambda x: x.split('|')[1:], fc_yy_rerun_str)

    print_message("re-running firstcal on XX pol", type=1)
    fc_xx_exits = M(fc_run, enumerate(fc_xx_rerun_str))
    print_message("re-running firstcal on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(fc_yy_rerun_str))

    fc_out.close()
    fc_err.close()
    print "fc_xx_exits:", fc_xx_exits
    print "fc_yy_exits:", fc_yy_exits

if rerun_fcmets is True:
    print_message("re-running firstcal metrics", type=1)

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
    fc_xx_rerun_ind = np.loadtxt(os.path.join(out_dir, 'fc_xx_rerun.tab'), dtype=np.str)
    fc_xx_rerun_ind = np.array(map(lambda x: x.split('|')[0], fc_xx_rerun_ind), np.int)
    fc_yy_rerun_ind = np.loadtxt(os.path.join(out_dir, 'fc_yy_rerun.tab'), dtype=np.str)
    fc_yy_rerun_ind = np.array(map(lambda x: x.split('|')[0], fc_yy_rerun_ind), np.int)

    print_message("running firstcal metrics on XX pol", type=1)
    fc_xx_exits = M(fc_run, enumerate(np.array(fc_xx_files)[fc_xx_rerun_ind]))
    print_message("running firstcal metrics on YY pol", type=1)
    fc_yy_exits = M(fc_run, enumerate(np.array(fc_yy_files)[fc_yy_rerun_ind]))

    fc_out.close()
    fc_err.close()
    print "fc_xx_exits:", fc_xx_exits
    print "fc_yy_exits:", fc_yy_exits

    # gather fc metrics
    xants_xx = np.loadtxt(os.path.join(out_dir, 'xants_xx.tab'), dtype=np.bool)
    xants_yy = np.loadtxt(os.path.join(out_dir, 'xants_yy.tab'), dtype=np.bool)
    fc_xx_xants = []
    fc_yy_xants = []
    fc_xx_rerun = []
    fc_yy_rerun = []

    for i, fcmxx in enumerate(np.array(fcm_xx_files)[fc_xx_rerun_ind]):
        xants = []
        if os.path.isfile(fcmxx) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fc_xx_files[i])
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmxx)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
            fc_xx_xants.append(xants)
            if len(xants) > 0:
                fc_xx_rerun.append(i)
        else:
            fc_xx_rerun.append(i)
            fc_xx_xants.append(xants)

    for i, fcmyy in enumerate(np.array(fcm_yy_files)[fc_yy_rerun_ind]):
        xants = []
        if os.path.isfile(fcmyy) is True:
            FC = hq.firstcal_metrics.FirstCal_Metrics(fc_yy_files[i])
            FC.plot_delays(save=True)
            plt.close()
            FC.load_metrics(fcmyy)
            ant_dly_std = np.array(FC.metrics['ant_std'].values())
            if ant_dly_std.max() > 0.5:
                xants.append(ants[np.argmax(ant_dly_std)])
            if len(FC.metrics['rot_ants']) > 0:
                xants.extend(FC.metrics['rot_ants'])
            if len(xants) > 0:
                fc_yy_rerun.append(i)
            fc_yy_xants.append(xants)
        else:
            fc_yy_rerun.append(i)
            fc_yy_xants.append(xants)

if collate_fcmets is True:
    print_message("collating firstcal metrics", type=1)

    mets_xx = OrderedDict({'file_times':[],
                            'ant_dly_std':OrderedDict(map(lambda a: (a, []), ants)),
                            'good_sol':[],
                            'rot_ants':[],
                            'xants':[]})
    mets_yy = OrderedDict({'file_times':[],
                            'ant_dly_std':OrderedDict(map(lambda a: (a, []), ants)),
                            'good_sol':[],
                            'rot_ants':[],
                            'xants':[]})

    for i, fcmxx in enumerate(fcm_xx_files):
        fc_xx_mets = hq.firstcal_metrics.load_firstcal_metrics(fcmxx)
        mets_xx['good_sol'].append(fc_xx_mets['good_sol'])
        mets_xx['rot_ants'].append(fc_xx_mets['rot_ants'])
        mets_xx['xants'].append(fc_xx_mets['bad_ants'])
        for a in ants:
            if a in fc_xx_mets['ants']:
                mets_xx['ant_dly_std'][a].append(fc_xx_mets['ant_std'][a])
            else:
                mets_xx['ant_dly_std'][a].append(np.nan)
    for i, fcmyy in enumerate(fcm_yy_files):
        fc_yy_mets = hq.firstcal_metrics.load_firstcal_metrics(fcmyy)
        mets_yy['good_sol'].append(fc_yy_mets['good_sol'])
        mets_yy['rot_ants'].append(fc_yy_mets['rot_ants'])
        mets_yy['xants'].append(fc_yy_mets['bad_ants'])
        for a in ants:
            if a in fc_yy_mets['ants']:
                mets_yy['ant_dly_std'][a].append(fc_yy_mets['ant_std'][a])
            else:
                mets_yy['ant_dly_std'][a].append(np.nan)

    # plot dly std
    fig = plt.figure(figsize=(14,8), dpi=100)
    ax = fig.add_subplot(111)
    ax.grid(True)
    # plotting
    cm_func = plt.get_cmap('nipy_spectral')
    cm = cm_func(np.linspace(0, 0.95, len(ants)))
    p = map(lambda x: ax.plot(np.array(jd_files, np.float), x[1], marker='o', ls='-', c=cm[x[0]]), enumerate(np.array(mets_xx['ant_dly_std'].values())))
    # axes
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    [t.set_rotation(20) for t in ax.get_xticklabels()]
    [t.set_fontsize(12) for t in ax.get_xticklabels()]
    [t.set_fontsize(12) for t in ax.get_yticklabels()]
    ax.set_xlabel('decimal of JD {}'.format(JD), fontsize=14)
    ax.set_ylabel('standard deviation [nanosec]', fontsize=16)
    ax.set_title("standard deviation of delay fluctuation per file for XX pol", fontsize=14)
    ax1 = fig.add_axes([.95,0.1,0.05,0.8])
    ax1.axis('off')
    leg = ax1.legend(np.concatenate(p), ants, ncol=1)
    fig.savefig('dlystd_xx_{}.png'.format(JD))
    plt.close()

    # plot dly std
    fig = plt.figure(figsize=(14,8), dpi=100)
    ax = fig.add_subplot(111)
    ax.grid(True)
    # plotting
    cm_func = plt.get_cmap('nipy_spectral')
    cm = cm_func(np.linspace(0, 0.95, len(ants)))
    p = map(lambda x: ax.plot(np.array(jd_files, np.float), x[1], marker='o', ls='-', c=cm[x[0]]), enumerate(np.array(mets_yy['ant_dly_std'].values())))
    # axes
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    [t.set_rotation(20) for t in ax.get_xticklabels()]
    [t.set_fontsize(12) for t in ax.get_xticklabels()]
    [t.set_fontsize(12) for t in ax.get_yticklabels()]
    ax.set_xlabel('decimal of JD {}'.format(JD), fontsize=14)
    ax.set_ylabel('standard deviation [nanosec]', fontsize=16)
    ax.set_title("standard deviation of delay fluctuation per file for YY pol", fontsize=14)
    ax1 = fig.add_axes([.95,0.1,0.05,0.8])
    ax1.axis('off')
    leg = ax1.legend(np.concatenate(p), ants, ncol=1)
    fig.savefig('dlystd_yy_{}.png'.format(JD))
    plt.close()


###################
### Run OmniCal ###
###################

# eliminate files
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
        pol = os.path.basename(uvfile).split('.')[3]
        cmd = "omni_run.py -C {} -p {} --ex_ants {} --overwrite --firstcal={} --omnipath={} {}".format(calfile, pol, xants, cffile, out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=oc_out, stderr=oc_err)
        return out

    # run oc
    print_message("running omnical on XX pol", type=1)
    oc_xx_exits = M(oc_run, enumerate(zip(fc_xx_files, xx_files, xants_xx)))
    print_message("running omnical on YY pol", type=1)
    oc_yy_exits = M(oc_run, enumerate(zip(fc_yy_files, yy_files, xants_yy)))

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
    oc_xx_exits = M(oc_run, enumerate(zip(fc_xx_files, oc_xx_files)))
    print_message("running omnical metrics on YY pol", type=1)
    oc_yy_exits = M(oc_run, enumerate(zip(fc_yy_files, oc_yy_files)))

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
    oa_xx_exits = M(oa_run, enumerate(zip(oc_xx_files, xx_files)))
    print_message("running omnical apply on YY pol", type=1)
    oa_yy_exits = M(oa_run, enumerate(zip(oc_yy_files, yy_files)))

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
        uvfile = inp[1]
        if i % (Nfiles//10) == 0: print "{:3d}%".format(int(np.round(float(i)/Nfiles, 1)*100))
        cmd = "xrfi_run.py --xrfi_path={} {}".format(out_dir, uvfile)
        out = subprocess.call(cmd, shell=True, stdout=rfi_out, stderr=rfi_err)
        return out

    # run ocm
    print_message("running xrfi  on XX pol", type=1)
    rfi_xx_exits = M(rfi_run, enumerate(xx_files))
    print_message("running xrfi on YY pol", type=1)
    rfi_yy_exits = M(rfi_run, enumerate(yy_files))

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
    rfi_xx_exits = M(rfi_run, enumerate(zip(flag_xx_files, oa_xx_files)))
    print_message("running rfi apply on YY pol", type=1)
    rfi_yy_exits = M(rfi_run, enumerate(zip(flag_yy_files, oa_yy_files)))

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


abs_xx_files = map(lambda x: x+"A", or_xx_files)
abs_yy_files = map(lambda x: x+"A", or_yy_files)
abs_data_path = "/data4/paper/SimulatedData/HERA19_GSM2008_NicCST2016_full_day/true_visibilities"

if abs_cal is True:
    print_message("absolute calibrating data", type=1)

    # get red bl groups
    for rbl in red_bls:
        if (20, 89) in rbl or (89, 20) in rbl:
            ew_bls = rbl
        elif (20, 81) in rbl or (81, 20) in rbl:
            ne_bls = rbl
        elif (20, 22) in rbl or (22, 20) in rbl:
            nw_bls = rbl
        elif (20, 112) in rbl or (112, 20) in rbl:
            ns_bls = rbl

    # get ephem object
    loc = aipy.cal.get_aa(calfile, np.array([0.15]))

    def jd2lst(loc, jd):
        loc.date = Time(jd, format='jd').datetime
        return loc.sidereal_time()

    # get abs data files
    sim_files = sorted(glob.glob(os.path.join(abs_data_path, "2*")))
    sim_xx_files = np.array(fnmatch.filter(sim_files, '*xx*'))
    sim_yy_files = np.array(fnmatch.filter(sim_files, '*yy*'))

    # get lst for each file
    sim_xx_lsts = np.array(map(lambda x: jd2lst(loc, float('.'.join(os.path.basename(x).split('.')[:2]))), sim_xx_files))
    sim_yy_lsts = np.array(map(lambda x: jd2lst(loc, float('.'.join(os.path.basename(x).split('.')[:2]))), sim_yy_files))

    # choose files to use for abs cal
    or_xx_lst_arr = np.array(map(lambda x: jd2lst(loc, float('.'.join(os.path.basename(x).split('.')[1:3]))), or_xx_files))*12/np.pi
    # fornax A at 2457680.5028028
    abscal_xx_files = np.array(sorted(np.array(or_xx_files)[np.argsort(np.abs(jd_files.astype(np.float)-2457680.5028028))][:2]))

    # iterate over files
    for i, xxfile in enumerate(abscal_xx_files):
        # load file data
        uvd = UVData()
        uvd.read_miriad(xxfile)
        data_freqs = uvd.freq_array.squeeze()
        data_lsts = np.unique(uvd.lst_array)

        # get file lst
        lst = jd2lst(loc, float('.'.join(xxfile.split('.')[1:3])))

        # get nearest two simulation files
        simfiles = sorted(sim_xx_files[np.argsort(np.abs(sim_xx_lsts-lst))][:2])

        # load sim data
        sim = UVData()
        sim.read_miriad(simfiles)
        sim_freqs = sim.freq_array.squeeze()
        sim_lsts = np.unique(sim.lst_array)

        # interpolate sim data onto file data coordinates
        sim_ew_data = interpolate.interp2d(sim_freqs, sim_lsts, np.real(sim.get_data(ew_bls[0])))(data_freqs, data_lsts).astype(np.complex128)
        sim_ew_data += 1j*interpolate.interp2d(sim_freqs, sim_lsts, np.imag(sim.get_data(ew_bls[0])))(data_freqs, data_lsts)
        sim_ew_amp = np.abs(sim_ew_data)
        sim_ew_phs = np.angle(sim_ew_data)
        sim_ne_data = interpolate.interp2d(sim_freqs, sim_lsts, np.real(sim.get_data(ne_bls[0])))(data_freqs, data_lsts).astype(np.complex128)
        sim_ne_data += 1j*interpolate.interp2d(sim_freqs, sim_lsts, np.imag(sim.get_data(ne_bls[0])))(data_freqs, data_lsts)
        sim_ne_amp = np.abs(sim_ne_data)
        sim_ne_phs = np.angle(sim_ne_data)
        sim_nw_data = interpolate.interp2d(sim_freqs, sim_lsts, np.real(sim.get_data(nw_bls[0])))(data_freqs, data_lsts).astype(np.complex128)
        sim_nw_data += 1j*interpolate.interp2d(sim_freqs, sim_lsts, np.imag(sim.get_data(nw_bls[0])))(data_freqs, data_lsts)
        sim_nw_amp = np.abs(sim_nw_data)
        sim_nw_phs = np.angle(sim_nw_data)
        sim_ns_data = interpolate.interp2d(sim_freqs, sim_lsts, np.real(sim.get_data(ns_bls[0])))(data_freqs, data_lsts).astype(np.complex128)
        sim_ns_data += 1j*interpolate.interp2d(sim_freqs, sim_lsts, np.imag(sim.get_data(ns_bls[0])))(data_freqs, data_lsts)
        sim_ns_amp = np.abs(sim_ns_data)
        sim_ns_phs = np.angle(sim_ns_data)
 
        # get badants
        xx_badants = np.unique((','.join(xants_xx)).split(',')).astype(int)
        yy_badants = np.unique((','.join(xants_yy)).split(',')).astype(int)

        # average red bls together
        ew_vis = np.median(np.array([uvd.get_data(bl) for bl in ew_bls if (bl[0] not in xx_badants) and (bl[1] not in xx_badants)]), axis=0)
        ew_amp = np.abs(ew_vis)
        ew_phs = np.angle(ew_vis)
        ne_vis = np.median(np.array([uvd.get_data(bl) for bl in ne_bls if (bl[0] not in xx_badants) and (bl[1] not in xx_badants)]), axis=0)
        ne_amp = np.abs(ne_vis)
        ne_phs = np.angle(ne_vis)
        nw_vis = np.median(np.array([uvd.get_data(bl) for bl in nw_bls if (bl[0] not in xx_badants) and (bl[1] not in xx_badants)]), axis=0)
        nw_amp = np.abs(nw_vis)
        nw_phs = np.angle(nw_vis)
        ns_vis = np.median(np.array([uvd.get_data(bl) for bl in ns_bls if (bl[0] not in xx_badants) and (bl[1] not in xx_badants)]), axis=0)
        ns_amp = np.abs(ns_vis)
        ns_phs = np.angle(ns_vis)

        # take ratios
        ew_ratio = sim_ew_data / ew_vis
        ne_ratio = sim_ne_data / ne_vis
        nw_ratio = sim_nw_data / nw_vis
        ns_ratio = sim_ns_data / ns_vis

        # make some plots
        def vis_plot(sim_vis, data_vis, bl='EW', fname=''):
            fig, axes = plt.subplots(2, 3, figsize=(18,10))
            axes = axes.ravel()
            fig.suptitle("{} baseline: {}".format(bl, fname))

            ax = axes[0]
            plt.sca(ax)
            uvt.plot.waterfall(sim_vis, mode='log', mx=1, drng=3)
            ax.set_title('simulation amplitude')
            plt.colorbar()
            ax = axes[1]
            plt.sca(ax)
            uvt.plot.waterfall(data_vis, mode='log', mx=1, drng=3)
            ax.set_title('data average amplitude')
            plt.colorbar()
            ax = axes[2]
            plt.sca(ax)
            uvt.plot.waterfall(sim_vis / data_vis, mode='log', mx=2, drng=3)
            ax.set_title('amplitude difference')
            plt.colorbar()
            ax = axes[3]
            plt.sca(ax)
            uvt.plot.waterfall(sim_vis, mode='phs', mx=np.pi, drng=2*np.pi)
            ax.set_title('simulation phase')
            plt.colorbar()
            ax = axes[4]
            plt.sca(ax)
            uvt.plot.waterfall(data_vis, mode='phs', mx=np.pi, drng=2*np.pi)
            ax.set_title('data average phase')
            plt.colorbar()
            ax = axes[5]
            plt.sca(ax)
            uvt.plot.waterfall(sim_vis / data_vis, mode='phs', mx=np.pi, drng=2*np.pi)
            ax.set_title('phase difference')
            plt.colorbar()

            plt.savefig('{}_vis.png'.format(bl))
            plt.close()

        vis_plot(sim_ew_data, ew_vis, bl='EW', fname=os.path.basename(xxfile))
        vis_plot(sim_ne_data, ne_vis, bl='NE', fname=os.path.basename(xxfile))
        vis_plot(sim_nw_data, nw_vis, bl='NW', fname=os.path.basename(xxfile))
        vis_plot(sim_ns_data, ns_vis, bl='NS', fname=os.path.basename(xxfile))


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













