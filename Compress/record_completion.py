#! /usr/bin/env python
"""
Record the completion of a task
1) create a new row in pdb.files
2) create a new row in pdb.history
3) delete the proper entry in pdb.proc
"""

from PDB import *
import sys
import optparse
import re

o = optparse.OptionParser()
o.add_option('-o', '--outfile', dest='outfile',type='string',default='',
        help='output file')
o.add_option('--host', dest='host', type='string',
        help='host of output file (for rsync) [Default `hostname`]')
o.add_option('-d', '--desc', dest='desc',type='string',
        help='Description of operation')
opts, args = o.parse_args()

def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]

infile = sys.argv[1]
outfile = opts.outfile

#update last_modified for infile
pdb.update('last_modified', "NOW()", 'files','filename',infile)

# Create a row for pdb.files
filecols = {}
filecols['JD'] = file2jd(infile)
filecols['filename'] = outfile

if not opts.host is None:
    filecols['host'] = opts.host
else:
    filecols['host'] = pdb.get('host', 'files', 'filename', infile)[0][0]

filecols['basefile'] = pdb.get('basefile', 'files', 'filename', infile)[0][0]
filecols['md5'] = __import__('hashlib').md5(outfile).hexdigest()
filecols['last_modified']="NOW()"

#reject if the md5 checksum doesn't agree.
insuffix  =  infile.split('.')[-1]
outsuffix = outfile.split('.')[-1]
if insuffix==outsuffix and not filecols['md5']==pdb.get('md5', 'files', 'filename', infile)[0][0]:
    print "md5 checksum failed. exiting"
    sys.exit(1)
pdb.addrow('files',filecols)

#add row to history table.
histcols = {}
histcols['input']=infile
histcols['output'] = outfile
histcols['operation'] = opts.desc
histcols['timestamp'] = "NOW()"
pdb.addrow('history',histcols)

#delete relevant entry of pdb.proc
pdb.delrow('proc',infile)
