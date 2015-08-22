#! /usr/bin/env python
"""
Add a row to pdb.hosts
syntax:
    add_host.py <hostname> <ipaddr> <username>

DFM
"""

from PDB import *
import sys

hostcols = {}
hostcols['hostname'], hostcols['IP'], hostcols['username'] = sys.argv[1:4]
if pdb.has_record('hosts', hostcols['hostname']):
    print "Record exists!"
else:
    pdb.addrow('hosts', hostcols)
