#! /usr/bin/env python
"""
Updates the outfile of history table to allow a competed record.
 -- exits 1 if bogus file is given.
"""

from PDB import *
import sys,optparse

o = optparse.OptionParser()
o.add_option('--log',type=str,
                        help='Log text. Usually saved in an env variable or catted from a file.')
opts, args = o.parse_args()
for outfile in args:
    if not pdb.has_record('history',outfile,col='output'):
        print "entry in the history column doesn't exist!"
        sys.exit(1)

    pdb.update('exit_status',"1",'history','output', outfile)
    pdb.update('stoptime',"NOW()",'history','output', outfile)
    pdb.update('log',opts.log,'history','output',outfile)
