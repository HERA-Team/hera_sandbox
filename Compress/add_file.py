#! /usr/bin/env python
"""
Add a row to pdb.files. Exits with 1 if a bogus host is used. Exits with 2 if the md5 checksum fails.
"""

from PDB import *
from socket import gethostname
import optparse
import re
import sys

o = optparse.OptionParser()
o.add_option('-i', '--infile', dest='infile',type='string',default='',
        help='input file')
opts, args = o.parse_args(sys.argv[1:])

def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]

hostname = gethostname()
infile  = opts.infile
outfile = args[0]

if not pdb.has_record('hosts', hostname):
    sys.exit(1)
if not pdb.has_record('files',infile):
    sys.exit(1)

filecols = {}
filecols['JD'] = file2jd(outfile)
filecols['filename'] = outfile
filecols['host'] = hostname
print pdb.get('basefile','files','filename',infile)
filecols['basefile'] = pdb.get('basefile','files','filename',infile)[0][0]
filecols['md5'] = __import__('hashlib').md5(outfile).hexdigest()
filecols['last_modified']="NOW()"

insuffix  =  infile.split('.')[-1]
outsuffix = outfile.split('.')[-1]
if insuffix==outsuffix and not filecols['md5']==pdb.get('md5', 'files','filename',infile)[0][0]:
    print "md5 checksum failed! exiting."
    sys.exit(2)

pdb.addrow('files', filecols)
