#! /usr/bin/env python
"""
Remove the row in the orders table with basename = basename of argv[1]

Be clever about it:
    --- leave orders whose staus is sys.argv[2] alone. (In normal usage, this will be 6-RSYNC and these will be read later to send data to USA
    --- if the order has already been deleted, write a new one with the status of sys.argv[2]. This should only happen for files on the edges of
    stills. In any case, maybe this should complain if that happens.

DFM
"""

from PDB import *
import sys

fname = sys.argv[1]
try:
    goodstatus = sys.argv[2]
except(IndexError):
    goodstatus = ''

basefile=pdb.get('basefile','files','filename',fname)[0][0]

#this is a really hacky workaround for having files on two different still machines.
try:
    status=pdb.get('status','orders','basefile',basefile)[0][0]
except(IndexError):
    shouldsend = False
    for suffix in 'DEF':
        shouldsend |= fname.endswith(suffix)
    if shouldsend:
        pdb.addrow('orders',{'status':'6-RSYNC','basefile':basefile})
        sys.exit(0)
    else:
        sys.exit(1)

if not status==goodstatus:
    pdb.delrow('orders', basefile)
