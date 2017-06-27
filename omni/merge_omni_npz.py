#!/usr/bin/env python
import numpy as np, capo.omni as omni, sys, os, optparse
"""
Combine omnical solutions for independent runs on xx and yy into
a single npz file per decimal JD. This lets us use omni_apply.py
to calibrate xx,xy,yx,yy at once, even though we only used info
from the linear instrumental polarizations.

Assumes xx and yy npz files are in the same directory.
Returns 
"""
o = optparse.OptionParser()
o.set_usage('merge_omni_npz.py [options] /path/to/*xx.npz')
opts,args = o.parse_args(sys.argv[1:])

def merge_dicts(*dict_args):
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

for xx_npz in args:
    yy_npz = xx_npz.replace("xx","yy")
    new_npz = xx_npz.replace("xx.","")
    print "Reading %s %s"%(xx_npz,yy_npz)    
    metax,gainsx,vismdlx,xtalkx = omni.from_npz(xx_npz)
    metay,gainsy,vismdly,xtalky = omni.from_npz(yy_npz)    
    meta = merge_dicts(metax,metay)
    gains = merge_dicts(gainsx,gainsy)
    vismdl = merge_dicts(vismdlx,vismdly)
    xtalk = merge_dicts(xtalkx,xtalky)
    print "    => %s"%new_npz
    omni.to_npz(new_npz,meta,gains,vismdl,xtalk)

