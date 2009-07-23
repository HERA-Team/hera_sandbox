#!/usr/bin/env python
#
#  init_ant_db.py
#  
#
#  Created by Danny Jacobs on 7/15/09.
#  PAPER Project
#
"""
A script for importing a list of antenna positions in a csv file into a 
capodb database.
usage
init_ant_db.py -o <output db filename> <input csv filename>

csv format:
Name,E,N,H
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
import capodb

o = optparse.OptionParser()
o.set_usage('init_ant_db.py [options] *.csv')
a.scripting.add_standard_options(o, cal=True)
o.add_option('-o','--dbfile', dest='filename', default='True',
    help='Filename of new database to be created.')
opts, args = o.parse_args(sys.argv[1:])
if len(args)==0: data = [s.split(',') for s in sys.stdin.readlines()]
else: data = [s.split(',') for s in open(args[0],'r').readlines()]
antpos = n.array([map(float,[d[1],d[2],d[3]]) for d in data])
antpos += -antpos[0]
print antpos
ants = []
beam = a.fit.Beam(n.linspace(0.1,0.2,num=256))
for i,pos in enumerate(antpos):
    ants.append(capodb.Antenna(pos[0],pos[1],pos[2],beam,id=i))
aa = capodb.AntennaArray([0,0],ants,filename=opts.filename) #location doesnt matter for db
aa.clear()
aa.save()