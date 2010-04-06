#! /usr/bin/env python

from matplotlib import use; use('agg')
import pylab as p
p.plot([0,2,1,3])
p.savefig('test.png')
