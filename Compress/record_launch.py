#! /usr/bin/env python
"""
Register a task in the proc table of pdb.
$ record_launch.py <full path to file> <task>
"""

from PDB import *
from socket import gethostname
import sys

hostname = gethostname()

if not pdb.has_record('hosts', hostname):
    sys.exit(1)

#update proc table
proccols = {}
proccols['host'] = hostname
proccols['filename'] = sys.argv[1]
proccols['operation'] = sys.argv[2]
proccols['starttime'] = "NOW()"
print proccols
pdb.addrow('proc', proccols)
