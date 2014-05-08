#! /usr/bin/env python
"""
Updates the outfile of history table to allow a competed record.
 -- exits 1 if bogus file is given.
"""

from PDB import *
import sys

outfile=sys.argv[1]

if not pdb.has_record('history',outfile,col='output'):
    print "entry in the history column doesn't exist!"
    sys.exit(1)

pdb.update('exit_status',"0",'history','output', outfile)
pdb.update('stoptime',"NOW()",'history','output', outfile)
