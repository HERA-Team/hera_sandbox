#! /usr/bin/env python

from PDB import *
import sys

thisfile = sys.argv[1]
basefile = pdb.get('basefile', 'files', 'filename', thisfile)[0][0]

provenance = []
while(thisfile != basefile):
    provenance.append(pdb.get('*','history','output', thisfile)[0])
    thisfile = provenance[-1][0]
provenance.append(pdb.get('*','history','output', thisfile)[0])

for infile,op,t,outfile in provenance[::-1]:
    print "%40s -> %10s -> %40s (%s)"%(infile, op, outfile, t)
