#! /usr/bin/env python

"""
Add a row to the pdb.observations table and populate the pdb.files tables accordingly.

DFM
"""

from initDB import pdb
import sys

hostname = None
list_of_jds = []

for arg in sys.argv[1:]:
    if arg.startswith('245'):
        list_of_jds.append(float(arg))
    else:
        hostname = arg

if not pdb['hosts'].has_record(hostname):
    sys.exit(1)

