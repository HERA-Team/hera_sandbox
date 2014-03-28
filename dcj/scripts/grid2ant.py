#! /usr/bin/env python
from capo.dfm import grid2ij
import optparse,sys
import aipy as a,numpy as  n
"""
input a grid spacing (i,j) and cal file
output a list of baselines

example:

"""

o = optparse.OptionParser()
o.set_usage('grid2int.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

#aa = a.cal.get_aa(opts.cal,n.array([0.15]))
exec("from %s import prms"%opts.cal)
out = ''
for arg in args:
    out += ','+grid2ij(prms['ant_layout'])[0][arg]
print out





