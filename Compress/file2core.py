#! /usr/bin/env python

import sys
import os

files=sys.argv[1:]
Nfiles=len(files)
Nproc=int(os.environ['Nproc'])
n=int(os.environ['n'])

if (Nproc <= Nfiles):
    NperProc=Nfiles/Nproc
    #print Nfiles,Nproc,NperProc
    print ' '.join(sorted(files[NperProc*n:NperProc*(n+1)]))
else:
    if n < Nfiles:
        print files[n]
    else:
        print ""
