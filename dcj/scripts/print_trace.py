#!/usr/bin/env python
#
#  print_trace.py
#  
#
#  Created by Danny Jacobs on 4/5/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

#o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
#opts, args = o.parse_args(sys.argv[1:])

def print_param_head(prms):
    keys = prms.keys()
    keys.sort()
    if type(prms[keys[0]]) is dict: keys = [keys[0]]; #assumes all records are the same.
    for k in keys:
        v = prms[k]
        if (type(v) is dict and v=={}) or v is None or \
            (type(v) is list and v == {}):
            continue
        if type(v)==dict:
            if 'jys' in prms[k].keys(): print 'src',
            elif 'pos' in prms[k].keys(): print 'ant',
            else: print k,    
            print_param_head(prms[k])
        else:
            print k,
    print
def print_param_line(prms):
    keys = prms.keys()
    keys.sort()
    for k in keys:
        v = prms[k]
        if (type(v) is dict and v=={}) or v is None or \
            (type(v) is list and v == {}):
            continue
        if type(v)==dict:
            print k,    
            print_param_line(prms[k])
        else:
            if not type(v) is list:
                try: v = [list(v)]
                except(TypeError): v = [v]
            print '\t'.join(map(str,v)),
    print


def print_params(prms, indent='', grad=None):
    """Print nice looking representation of a parameter dictionary."""
    keys = prms.keys()
    keys.sort()
    for k in keys:
        v = prms[k]
        if (type(v) is dict and v == {}) or v is None or \
                (type(v) is list and v == []):
            continue
        if type(v) == dict:
            print indent, k
            if grad is None: print_params(v, indent + '  ')
            else: print_params(v, indent + '  ', grad[k])
        else:
            print indent, k,
            if grad is None:
                if not type(v) is list:
                    try: v = [list(v)]
                    except(TypeError): v = [v]
                if len(v) == 1: print v[0]
                else:
                    print
                    for i in v: print indent, ' ', i
            else:
                print
                print indent, v, '\t<', grad[k], '>'
                if not type(v) is list:
                    try: v = [list(v)]
                    except(TypeError): v = [v]
                for i in len(v):
                    print indent, ' ', v[i], '\t<', grad[k][i], '>'
