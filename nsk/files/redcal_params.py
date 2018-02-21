"""
redcal_params.py
----------------

parameter file for redcal_pipeline.py
"""
import os
import sys
import glob
import pathos

## set flags ##
T = True
F = False

overwrite           = T

run_firstcal        = F
run_fcmets          = F
collate_fcmets      = F

run_omnical         = T
run_ocmets          = T
apply_omni          = T

run_rfi             = T
apply_rfi           = T
plot_reds           = F

multiprocess        = F
Nproc               = 16

def print_message(message, type=0):
    if type == 0:
        print message
    elif type == 1:
        print "\n"+message+"\n"+"-"*40

# get JD
JD = 2458042

# assign mp pooling
if multiprocess is True:
    pool = pathos.multiprocessing.Pool(Nproc)
    M = pool.map
else:
    M = map

# Get data files
data_path = os.path.join("/lustre/aoc/projects/hera/H1C_IDR1", str(JD))
out_dir = os.path.join("/lustre/aoc/projects/hera/nkern/data/H1C", str(JD))

