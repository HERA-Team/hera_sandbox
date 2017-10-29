# -*- mode: python; coding: utf-8 -*-
# Copyright (c) 2017 The HERA Collaboration

from __future__ import print_function, division, absolute_import
import os
import re
import sys
import time
import ConfigParser as configparser
from configparser import ConfigParser, ExtendedInterpolation


def get_config_entry(config, header, item, required=True):
    '''
    Helper function to extract specific entry from config file.
    '''
    if config.has_option(header, item):
        entries = config.get(header, item).split(',')
        entries = [entry.strip() for entry in entries]
        return entries
    else:
        if not required:
            return []
        else:
            raise AttributeError("Error processing config file: item \"{0}\" under header \"{1}\" is "
                                 "required, but not specified".format(item, header))


def make_outfile_name(obsid, action, pol_list=[]):
    '''
    Make a list of unique output files names for each stage and polarization.
    '''
    outfiles = []
    if len(pol_list) > 1:
        for pol in pol_list:
            of = "{0}.{1}.{2}.out".format(obsid, action, pol)
            outfiles.append(of)
    else:
        of = "{0}.{1}.out".format(obsid, action)
        outfiles.append(of)
    return outfiles


def prep_args(args, obsid, pol):
    '''
    Substitute the polarization string in a filename/obsid with the specified one.
    '''
    # replace pol if present
    match = re.search(r'zen\.\d{7}\.\d{5}\.(.+)\.', obsid)
    if match:
        obs_pol = match.group(1)
        basename = re.sub(obs_pol, pol, obsid)
    else:
        basename = obsid
    # replace {basename} with actual basename
    return re.sub(r'\{basename\}', basename, args)


def build_makeflow_from_config(obsids, config_file, mf_name=None):
    '''
    Construct a makeflow file from a config file.

    Config file structure:

    [STAGENAME]
    prereqs = STAGENAME1, STAGENAME2
    args = arg1, arg2
    ncpu = 1
    mem = 5000 (MB)


    Mini-language for replacement:
    basename = "zen.2458000.12345.uv"

    '''
    # read in config file
    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(config_file) 
    workflow = get_config_entry(config, 'WorkFlow', 'actions')

    # get general options
    pol_list = get_config_entry(config, 'Options', 'pols', required=False)
    if len(pol_list) == 0:
        # make a dummy list of length 1, to ensure we perform actions later
        pol_list = ['']
    path_to_do_scripts = get_config_entry(config, 'Options', 'path_to_do_scripts')[0]
    parent_dirs = get_config_entry(config, 'Options', 'parent_dirs')[0]

    # open file for writing
    cf = os.path.basename(config_file)
    if mf_name is not None:
        fn = "{0}.{1}.mf".format(mf_name, cf)
    else:
        t = str(int(time.time()))
        fn = "{0}.{1}.mf".format(t, cf)

    # write file
    with open(fn, "w") as f:
        for obsid in obsids:
            for ia, action in enumerate(workflow):
                # start list of input files
                infiles = []

                # get dependencies
                prereqs = get_config_entry(config, action, "prereqs", required=False)
                if len(prereqs) > 0:
                    for prereq in prereqs:
                        try:
                            ip = workflow.index(prereq)
                        except ValueError:
                            raise ValueError("Prereq \"{0}\" for action \"{1}\" not found in main "
                                             "workflow".format(prereq, action))
                        outfiles = make_outfile_name(obsid, prereq, pol_list)
                        for of in outfiles:
                            infiles.append(of)

                # add command to infile list
                # this implicitly checks that do_{STAGENAME}.sh script exists
                command = "do_{}.sh".format(action)
                command = os.path.join(path_to_do_scripts, command)
                infiles.append(command)

                # also add previous outfiles to input requirements
                if ia > 0:
                    for of in outfiles_prev:
                        infiles.append(of)
                infiles = ' '.join(infiles)

                # make argument list
                args = get_config_entry(config, action, "args", required=False)
                args = ' '.join(args)

                # make outfile name
                outfiles = make_outfile_name(obsid, action, pol_list)

                # get directory where makeflow file is (since makeflow looks here for output files)
                cwd = os.getcwd()

                # make rules
                for pol, outfile in zip(pol_list, outfiles):
                    # replace '{basename}' with actual filename
                    # replace polarization string
                    prepped_args = prep_args(args, obsid, pol)

                    # make logfile name
                    # logfile will capture stdout and stderr
                    logfile = re.sub('\.out', '.log', outfile)

                    # make a small wrapper script that will run the actual command
                    # can't embed if; then statements in makeflow script
                    wrapper_script = re.sub('\.out', '.sh', outfile)
                    wrapper_script = "wrapper_{}".format(wrapper_script)
                    with open(wrapper_script, "w") as f2:
                        print("#!/bin/bash", file=f2)
                        print("source activate hera", file=f2)
                        print("date", file=f2)
                        print("cd {}".format(parent_dirs), file=f2)
                        print("{0} {1}".format(command, prepped_args), file=f2)
                        print("if [ $? -eq 0 ]; then", file=f2)
                        print("  cd {}".format(cwd), file=f2)
                        print("  touch {}".format(outfile), file=f2)
                        print("fi", file=f2)
                        print("date", file=f2)
                    # make file executable
                    os.chmod(wrapper_script, 0o755)

                    # first line lists target file to make (dummy output file), and requirements
                    # second line is "build rule", which runs the shell script and makes the output file
                    line1 = "{0}: {1}".format(outfile, infiles)
                    line2 = "\t./{0} > {1} 2>&1\n".format(wrapper_script, logfile)
                    print(line1, file=f)
                    print(line2, file=f)

                # save previous outfiles for next time
                outfiles_prev = outfiles

    return


