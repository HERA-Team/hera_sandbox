#! /usr/bin/env python

import optparse
from numpy import pi, arange
import math
import sys
import os


o = optparse.OptionParser()
o.add_option('--lst_rng', dest='rng', default='0_24',
    help="Total range of LST's to bin (hours) Default=0_24")
o.add_option('--file_len', dest='flen', type=float, default=0.5,
    help='Lenght of a file in hours of LST. Default=0.5')
opts, args = o.parse_args(sys.argv[1:])

#figure out which files I'm going to eventually make
lststart,lststop = map(float, opts.rng.split('_'))
lsts = arange(lststart,lststop + opts.flen,opts.flen)
lsts *= pi/12. #convert hours to radians.
Nfiles = len(lsts) - 1
lstfiles = []
for i in range(Nfiles):
    lstfiles.append("%f_%f"%(lsts[i],lsts[i+1]))

#where do i fit in the gridengine environment?
task1    = int(os.environ['SGE_TASK_FIRST'])
taskn    = int(os.environ['SGE_TASK_LAST'])
i = int(os.environ['SGE_TASK_ID']) - 1
Ntasks = taskn-task1

if Ntasks <= Nfiles:
    num = int(math.ceil(float(Nfiles)) / (Ntasks + 1))
    print ' '.join(lstfiles[num*i:num*(i+1)])
else:
    print lstfiles[i%Nfiles]
