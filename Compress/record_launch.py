#! /usr/bin/env python
"""
Register a task in the proc table of pdb.
$ record_launch.py <full path to file> <task>
"""

from PDB import *
from socket import gethostname
import optparse
import sys

o = optparse.OptionParser()
o.add_option('-i', '--infile', dest='infile',type='string',default='',
        help='input file')
o.add_option('-d','--desc',dest='desc', type='string',
        help='Description of operation')
opts, args = o.parse_args()

hostname = gethostname()
infile = opts.infile
outfile = args[0]

if not pdb.has_record('hosts', hostname):
    sys.exit(1)
if not pdb.has_record('files',infile):
    sys.exit(1)

pdb.update('last_modified',"NOW()",'files','filename',infile)

#update history table
histcols = {}
histcols['input']  = infile
histcols['output'] = outfile
histcols['operation'] = opts.desc
histcols['starttime'] = "NOW()"
pdb.addrow('history', histcols)
