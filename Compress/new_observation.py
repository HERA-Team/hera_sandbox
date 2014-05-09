#! /usr/bin/env python

"""
Add a row to the pdb.observations table and populate the pdb.files tables accordingly.
Add a row to pdb.orders indicating that compression needs done.

DFM
"""

from PDB import *
from socket import gethostname
import sys
import re

hostname = gethostname()
prefix = None
list_of_jds = []
new_files=["%s:%s"%(hostname,a) for a in sys.argv[1:]]

def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]

for arg in sys.argv[1:]:
    if 'zen.245' in arg:
        jd = file2jd(arg)
        if prefix is None:
            prefix = '/'.join(arg.split('/')[:-1])
        if not jd in list_of_jds:
            list_of_jds.append(jd)

list_of_jds.sort()

#exit if the host is bad.
if not pdb.has_record('hosts',hostname):
    print "No entry for host %s exists. Exiting (1)"%hostname
    sys.exit(1)

for i,jd in enumerate(list_of_jds):
    obscols = {}
    obscols['JD'] = jd
    for pi in 'xy':
        for pj in 'xy':
            pol_fname = "%s:%s/zen.%s.%s.uv"%(hostname,prefix,jd,pi+pj)
            if pol_fname in new_files:
                obscols[pi+pj] = pol_fname
                # add an entry to pdb.files here. Otherwise observations will break.
                filecols = {}
                filecols['JD'] = jd
                filecols['filename'] = pol_fname
                filecols['basefile'] = pol_fname
                filecols['host'] = hostname
                filecols['created_on'] = "NOW()"
                filecols['last_modified'] = "NOW()"
                pdb.addrow('files',filecols)
                # update history to include created file.
                histcols = {}
                histcols['input'] = filecols['filename']
                histcols['output'] = filecols['filename']
                histcols['operation'] = "0-CREATED"
                histcols['host'] = hostname
                histcols['starttime'] = "NOW()"
                histcols['stoptime'] = "NOW()"
                histcols['exit_status'] = '0'
                pdb.addrow('history',histcols)

                ordercols = {}
                ordercols['filename'] = pol_fname
                ordercols['status'] = 'CREATED'
                pdb.addrow('orders', ordercols)
    try:
        obscols['jd_hi'] = list_of_jds[i+1]
    except(IndexError):
        pass
    if i >= 1:
        obscols['jd_lo'] = list_of_jds[i-1]

    obscols['created_on'] = "NOW()"
    obscols['last_modified'] = "NOW()"

    pdb.addrow('observations',obscols)
