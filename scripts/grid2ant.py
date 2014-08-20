#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from capo.dfm import grid2ij
import optparse,sys
import aipy as a,numpy as  n
"""
input a grid spacing (i,j) and cal file
output a list of baselines

example:
grid2ant.py -C psa898_v003 --seps="0,1;-1,1"

"""

o = optparse.OptionParser()
o.set_usage('grid2int.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--seps',type=str,
   help="The seperation types you want. ; delimited")
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

#aa = a.cal.get_aa(opts.cal,n.array([0.15]))
exec("from %s import prms"%opts.cal)
out = ''
for arg in opts.seps.split(';'):
    out += ','+grid2ij(prms['ant_layout'])[0][arg]
print out[1:]





