#! /usr/bin/env python
"""
Register a task in the proc table of pdb.
$ record_launch.py <full path to file> <task>
"""

from PDB import *
import optparse
import sys

o = optparse.OptionParser()
o.add_option('-i', '--infile', dest='infile',type='string',default='',
        help='input file')
o.add_option('-d','--desc',dest='desc', type='string',
        help='Description of operation')
opts, args = o.parse_args()

hostname = gethostname()
infile  = opts.infile
outfile = args[0]

if not pdb.has_record('hosts', hostname):
    if pdb.verbose:
        print "Unidentified host %s, exiting (1)"%hostname
    sys.exit(1)
if not pdb.has_record('files',infile):
    if pdb.verbose:
        print "Unidentified file %s, exiting (1)"%filename
    sys.exit(1)
pdb.update('last_modified',"NOW()",'files','filename',infile)

#update history table
histcols = {}
histcols['input']  = infile
histcols['output'] = outfile
histcols['host'] = hostname
histcols['operation'] = opts.desc
histcols['starttime'] = "NOW()"
#get the infile basefile, note that we already checked that it exists in the files table.
histcols['basefile'] = pdb.get('basefile','files','filename',infile)[0][0] #get basefile from files where filename=infile

#histcols['log'] = whatever
pdb.addrow('history', histcols)

