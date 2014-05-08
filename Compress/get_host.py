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
stdout=pdb.get('host','files','filename',fname)[0][0]
print stdout
