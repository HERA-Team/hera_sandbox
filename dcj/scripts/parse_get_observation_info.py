#! /usr/bin/env python
import sys
lines = open(sys.argv[-1]).readlines()
newrec=0
for line in lines:
    if line.startswith('PicA_141'):
        newrec=1
        filename = line.split()[2]

    if newrec and line.startswith('(Az'):
        El = float(line.split('=')[1].split(',')[1].split(')')[0].strip())
        if El>60:
            print filename,El
        newrec=0
