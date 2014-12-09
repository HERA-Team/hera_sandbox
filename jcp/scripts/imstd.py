#! /usr/bin/env python

import aipy as a, numpy as n, sys

for file in sys.argv[1:]:
    im,kwd = a.img.from_fits(file)
    print file, n.round(n.std(im[400:600,400:600])*1000), 'mjy'
