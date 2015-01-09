#! /usr/bin/env python

import numpy as np, FortranFile, sys

inFname = sys.argv[1]
outFname = sys.argv[2]
f = FortranFile.FortranFile(inFname)
vector = f.readReals("d")
np.save(outFname,vector)
