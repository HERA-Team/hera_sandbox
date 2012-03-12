#!/usr/bin/env python
#
#  cat_stats.py
#  
#
#  Created by Danny Jacobs on 3/19/10.
#  PAPER Project
#
"""
Count matches between catalogs, disjoint sets, etc.
"""
import aipy as a, numpy as n,math as m,warnings
import sys, optparse,ephem,vo.table as vot, os,logging
from pylab import *
from astrogrid import ConeSearch
from cStringIO import StringIO


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)

opts, args = o.parse_args(sys.argv[1:])

def find_src_in_cat(src,cat,sep):
    """
    Given a src (src = aipy.amp.RadioFixedObject or lower) 
    and a catalog (cat = aipy.amp.RadioCatalog or lower), find the closest
    source in cat that falls within sep [deg].
    Return: closest source object and seperation (ephem.Angle).
    If no source found within sep, return None,closest_sep
    """
    last_sep=n.pi
    closest_src=None
    for nk,test_src in cat.iteritems():
        cursep = ephem.separation(src, test_src)
        if cursep <= last_sep:
            last_sep = cursep
            closest_src = test_src
    if last_sep<sep * a.img.deg2rad: return closest_src,last_sep
    else: return None,last_sep
def union_cats(catA,catB,sep,detailed=False):
    """
    Given two catalogs of sources, return a set of paired sources that are
    seperated by sep degrees or less. Detailed also returns vector of 
    seperations.
    union = union_cats(catA,catB,sep)
    union,seps = union_cats(catA,catB,sep,detailed=True)
     """
    print len(catA),len(catB)
    union = []
    seps = []
    for k in catA.keys():
        catA_src = catA[k]
        catB_src,last_sep = find_src_in_cat(catA_src,catB,sep)
#        if k=='060633-202201': print catB_src
        if last_sep<sep * a.img.deg2rad:
            union.append([catA_src,catB_src])
            seps.append(last_sep)
    if detailed: return union,seps
    else: return union


#setup logging
if opts.vverb: logging.basicConfig(level=logging.DEBUG)
elif opts.verb: logging.basicConfig(level=logging.INFO)
else: logging.basicConfig(level=logging.WARNING)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)    
    
    