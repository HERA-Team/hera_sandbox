#! /usr/bin/env python
"""
input an index i, a number of nodes, and a number of files to distribute
return a list of indices into any list of nodes

This is a "multicast" algorithm. If a file is sent to the edge slot on a host, it is also sent to the next host in the list
"""


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
