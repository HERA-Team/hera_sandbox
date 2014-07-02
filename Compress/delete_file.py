#! /usr/bin/env python
"""
Remove a row from pdb.files.
"""

from PDB import *
from socket import gethostname
import sys

delfile = sys.argv[1]

if not pdb.has_record('files', delfile):
    print 'No Record, exiting.'
    sys.exit(1)
else:
   pdb.delrow('files', delfile)

