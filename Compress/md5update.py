#! /usr/bin/env python
"""
Compute the md5 checksum and populate the md5 column of the pdb.files table.
"""

from PDB import *
import sys
import hashlib

for filename in sys.argv[1:]:
    md5 = hashlib.md5(filename).hexdigest()
    q = "UPDATE files SET md5='%s' WHERE filename='%s';"%(md5,filename)
    pdb.db.query(q)
