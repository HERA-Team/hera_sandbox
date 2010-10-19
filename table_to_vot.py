#!/usr/bin/env python
#
#  table_to_vot.py
#  
#
#  Created by Danny Jacobs on 3/17/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,atpy

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--add',type=str,
    help="Add variables to be applied across the board. ie --add mfreq:0.15,seq:100")
opts, args = o.parse_args(sys.argv[1:])


for file in args:
    lines = open(file).readlines()
    t = atpy.Table()
    t.add_column('Name',[s.split('\t')[0].strip() for s in lines if not s.startswith('#')],dtype=n.str)
    t.add_column('Ra',[s.split('\t')[1].strip() for s in lines if not s.startswith('#')],dtype=n.str)
    t.add_column('Dec',[s.split('\t')[2].strip() for s in lines if not s.startswith('#')],dtype=n.str)
    t.add_column('S_nu_',[s.split('\t')[3].strip() for s in lines if not s.startswith('#')],dtype='<f8')
    try: t.add_column('e_S_nu_',[s.split('\t')[4].strip() for s in lines if not s.startswith('#')],dtype='<f8')
    except(IndexError): t.add_column('e_S_nu_',['0' for s in lines if not s.startswith('#')],dtype='<f8')
    if not opts.add is None:
        adds = opts.add.split(',')
        for add in adds:
            t.add_column(add.split(':')[0],[add.split(':')[1] for s in lines if not s.startswith('#')],dtype='S14')
    print file + ' > ' + '.'.join(file.split('.')[:-1])+'.vot'
    t.write('.'.join(file.split('.')[:-1])+'.vot')
    