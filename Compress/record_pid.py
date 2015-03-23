#! /usr/bin/env python
"""
Register a task in the proc table of pdb.
$ record_pid.py <full path to file> --pid=<pid>
"""

from PDB import *
import optparse
import sys

o = optparse.OptionParser()
o.add_option('--pid',type=int,
                help='pid of running process')
opts, args = o.parse_args()
output=sys.argv[1]

if opts.pid is None:
    print "NO PID ENTERED. EXITING"
    sys.exit(1)
pdb.update('pid',opts.pid,'history','output',output)
