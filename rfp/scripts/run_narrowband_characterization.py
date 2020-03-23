# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 the HERA Collaboration
# Licensed under the 2-clause BSD license.

import glob
import multiprocessing as mp
import os
import re
import warnings

from pyuvdata import UVData
from rids.features import narrowband

# overview

# figure out IO details
# extract analysis parameters
### easiest solution: use yaml files to pull the above info TODO
### probably good practice: make an argparser to let the user TODO
### use this script if they aren't familiar with making yaml files
### can configure the module argparser to accept a config file

# make new UVData objects that contain only desired info
# (e.g. only autos if not using cross-correlations)
### make a module in narrowband to handle this XXX DONE
### configure script to parallelize this task using multiprocessing

# run characterization routine, parallelizing over files where appropriate
### again, use a module in narrowband for doing this XXX DONE
### parallelization currently must be handled at the script-level
### it is also unclear if parallelization over files is the most 
### efficient way to handle this, since xrfi has some native multiprocessing
### that speeds up part of the rfi station extraction step

# save the desired data products
# write out all the info for the rids file
# save the rids file
### make sure to subclass Rids to handle this, again in narrowband TODO
### write a function (maybe in data_management) that takes the output from 
### the processing pipeline and creates a rids file TODO

# exit
