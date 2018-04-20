#!/usr/bin/env python2.7
'''Script for checking the status of a HERA makeflow pipeline.
By Josh Dillon, jsdillon@berkeley.edu.'''

import numpy as np
import glob
from hera_opm.mf_tools import get_config_entry
import ConfigParser as configparser
from configparser import ConfigParser, ExtendedInterpolation
import os
from datetime import datetime
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Check the status of a pipeline. Prints out the total number of jobs of each " +
                                 "task in the workflow (by the wrapper*.sh files), the number completed (by the .out " +
                                 "files), the number currently running (which have starting but not stopping times in the" +
                                 ".log files), and the number errored (which have stopping times in the log but no .out" +
                                 "file). Also prints the average time elapsed for non-trivial jobs (i.e. those that take" + 
                                 "more than a second).")
parser.add_argument("--config_file", type=str, required=True, 
                    help="Absolute path to makeflow .cfg file.")
parser.add_argument("--working_dir", nargs='*', type=str, required=True,
                    help="Absolute path to pipeline working directory (or directories using *).")
args = parser.parse_args()
if np.all([not os.path.isdir(wdir) for wdir in args.working_dir]):
    raise ValueError('You must supply at least one directory using --working_dir')

# Read makeflow config file     
config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(args.config_file)
workflow = get_config_entry(config, 'WorkFlow', 'actions')


def elapsed_time(log_lines):
    '''Take a list of lines in a log file and calculates the elapsed time in minutes.
    
    Arguments:
        log_lines: list of lines in a log file, like that produced by [file].readlines()
    
    Returns:
        runtime: runtime in minutes. 
            If the file has a start time but no end time, returns -1
            If the file has no start time, returns -2
    '''
    try:
        start = datetime.strptime(log_lines[0], '%a %b %d %H:%M:%S %Z %Y ')
    except:
        start = None
    try:
        end = datetime.strptime(log_lines[-1], '%a %b %d %H:%M:%S %Z %Y ')
    except:
        end = None
    
    if (start is not None) and (end is None):
        return -1
    elif (start is None):
        return -2
    else:
        return ((end - start).seconds + 24.0*60*60*(end - start).days)/60.0

    
def inspect_log_files(log_files, out_files):
    '''Look at log files, compute the average non-zero runtime, and print example errors.
    
    Arguments:
        log_files: list of .log files
        out_files: list of .out files
    
    Returns:
        average_runtime: average runtime of all finished jobs that took longer than a second, in minutes
        nRunning: number of jobs believed to be running (start time with no stop time)
    '''
    error_warned = False
    runtimes = []
    for log_file in log_files:
        with open(log_file,'r') as f:
            log_lines = f.readlines()
        runtimes.append(elapsed_time(log_lines))
        
        #If it completed running but there's no .out file, it's suspected to have errored
        if (runtimes[-1] > 0) and (log_file.replace('.log','.out') not in out_files) and not error_warned:
            print '\n\nError Suspected (no .out found) in', log_file
            print '------------------------------------------------\n'
            print '\n'.join(log_lines)
            print '------------------------------------------------\n\n'
            error_warned = True
        
    runtimes = np.array(runtimes)
    finished_runtimes = runtimes[runtimes>1/60.0]    
    if len(finished_runtimes) > 0:
        average_runtime = np.mean(finished_runtimes)
    else:
        average_runtime = np.nan
       
    return average_runtime, np.sum(runtimes==-1)


# Run pipeline report
print '---------------\nPIPELINE REPORT\n---------------\n'
for job in workflow:
    print job + ':'
    total, logged, done = [], [], []
    for wdir in args.working_dir:
        if os.path.isdir(wdir):
            total += glob.glob(os.path.join(wdir,'wrapper_*.' + job + '.*sh'))
            logged += glob.glob(os.path.join(wdir,'*.' + job + '.*log'))
            done += glob.glob(os.path.join(wdir,'*.' + job + '.*out'))
    average_runtime, nRunning = inspect_log_files(logged, done)
    nErrored = len(logged) - (len(done) + nRunning)
    
    print 'Average (non-zero) runtime:', average_runtime, 'minutes' 
    print len(total), '\t|\t', len(done), '\t|\t', nRunning, '\t|\t', nErrored
    print 'total\t|\tdone\t|\trunning\t|\terrored\n\n'
    