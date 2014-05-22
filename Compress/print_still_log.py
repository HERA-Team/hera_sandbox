#! /usr/bin/env python
"""
Print the operation log for the input file.
Default operation is to print only the log for the most recent step.
--all prints all logs leading up to this file as well
usage:
    print_still_log.py host:filepath
examples:
    print_still_log.py pot0:/data0/jd2456566/zen.2456566.uvcRRE #prints the most recent thing to happen to this file. Probably the move from the still
    print_still_log.py still4:/data/zen.2456566.uvcRE --all #prints every log we got on this observation leading up to this file
"""
from PDB import *
import optparse,sys,os

o = optparse.OptionParser()
o.add_option('--all',action='store_true',
                        help='Print ')
opts, args = o.parse_args()

#basefile=pdb.get('basefile','files','filename',fname)[0][0]
#destination_dir = pdb.get('output','history',['operation','basefile'],['1-RSYNC',basefile])[0][0]
for filename in args:
    if opts.all:
        try:
            basefile = pdb.get('basefile','history','output',filename)[0][0]
        except(IndexError):
            print filename+"not found in database"
            continue
        log = pdb.get('log','history','basefile',basefile)
        ops = pdb.get('operation','history','basefile',basefile)
        for i,l in enumerate(log):
            if not l is None:
                print "==="*5
                print ops[i]
                print l[0]

    else:
        log = pdb.get('log','history','output',filename)[0][0]
        print log

