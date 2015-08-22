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
outfile = args[0]
if not pdb.has_record('history',outfile,col='output'):
    print "entry in the history column doesn't exist!"
    sys.exit(1)

pdb.update('exit_status',"0",'history','output', outfile)
pdb.update('stoptime',"NOW()",'history','output', outfile)
if not opts.log is None:
    pdb.update('log',opts.log,'history','output',outfile)
#update orders.
basefile=pdb.get('basefile','files','filename',outfile)[0][0]
op=pdb.get('operation','history','output',outfile)[0][0]
pdb.update('status',op,'orders','basefile',basefile)
