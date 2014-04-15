#! /usr/bin/env python
"""
Add a row to pdb.hosts
syntax:
    add_host.py <hostname> <ipaddr> <username>

DFM
"""

from initDB import pdb
import sys

cols = {}
cols['hostname'], cols['IP'], cols['username'] = sys.argv[1:4]
if pdb.has_record('hosts', cols['hostname']):
    print "Record exists!"
else:
    pdb.addrow('hosts', cols)
