#! /usr/bin/env python

import sys

i = int(sys.argv[1])
Nnodes = int(sys.argv[2])
Nfiles = int(sys.argv[3])

Nfiles_per_node=Nfiles/Nnodes
stillno=i/Nfiles_per_node
stillind=i%Nfiles_per_node

stdout=[stillno]

if (stillind == 0) and (stillno != 0):
    stdout += [stillno-1]

if (stillind + 1 == Nfiles_per_node) and (stillno +1 != Nnodes):
    stdout += [stillno+1]

print ' '.join(sorted(map(str,stdout)))
