#! /usr/bin/env python
"""
Pass a filename into sys.argv[1], return the host where that files lives.
   --- known problem: what if different files share the same filename?

DFM
"""

from PDB import *
import sys

pdb.verbose=False #just in case...
fname=sys.argv[1]
basefile=pdb.get('basefile','files','filename',fname)[0][0]
destination_dir = pdb.get('output','history',['operation','basefile'],['1-RSYNC',basefile])[0][0]

print destination_dir.split(':')[0]
