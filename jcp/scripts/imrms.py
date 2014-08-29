#! /usr/bin/env python

import aipy as a, numpy as n, sys

sum, cnt = 0,0
for file in sys.argv[1:]:
    print file
    im,kwd = a.img.from_fits(file)
    sum+=im**2
    cnt+=1

im=n.sqrt(sum/cnt)
a.img.to_fits('test.fits',im,**kwd)
