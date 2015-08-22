#! /usr/bin/env python
"""
Compute the md5 checksum and populate the md5 column of the pdb.files table.
"""

from PDB import *
import sys

for filename in sys.argv[1:]:
    md5 = md5sum(filename)
    q = "UPDATE files SET md5='%s' WHERE filename='%s';"%(md5,filename)
    if pdb.verbose: print q
    pdb.db.query(q)
