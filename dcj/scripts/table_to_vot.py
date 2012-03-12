#!/usr/bin/env python
#
#  table_to_vot.py
#  
#
#  Created by Danny Jacobs on 3/17/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,atpy,ephem as e

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--add',type=str,
    help="Add variables to be applied across the board. ie --add mfreq:0.15,seq:100")
o.add_option('-d',type=str,default='\t',
    help="Delimeter [default=\\t]")
opts, args = o.parse_args(sys.argv[1:])
if opts.d is '': opts.d = None
print 'Seperator between colons :%s:'%opts.d
types = {'ra':str,'dec':str,'name':str,'npix':int,'_ra':float,'_dec':float,
's_nu_':float,'s_nu_int_':float,'dr':float}
ucds = {"pos.eq.ra;meta.main":'_Ra',
        "pos.eq.dec;meta.main":'_Dec',
        "phot.flux.density":'S_nu_',
        "stat.error":'e_S_nu_',
        "meta.id;meta.main":'Name',
        "meta.record":'Seq',
        "spect.index":'n',
        "em.freq":'nu'}
ucds = dict(zip(ucds.values(),ucds.keys()))
for file in args:
    lines = open(file).readlines()
    lines = [l for l in lines if not l.startswith('#')]
    names = [s.strip() for s in lines[0].split(opts.d)]
    lines = lines[1:]
    print "columns:",names
    t = atpy.Table()
    print [len(s.split(opts.d)) for s in lines]
    for i,name in enumerate(names):
        try:
            col =  [s.split(opts.d)[i].strip() for s in lines]
#        col = [s.split('\t')[i].strip() for s in lines]
        except(IndexError):
            #Column is "irregular" flag column, 
            #this kinda sucks and only works for the last column
            col = []
            for l in lines:
                l = l.split()
                if len(names)== len(l):
                    col.append(l[i].strip())
                else: col.append('')
        try: 
            float(col[0])
            tipe = '<f8'
            if name.endswith('px'): type='i'
        except(ValueError):tipe = n.str
        if name in ucds.keys():
            ucd = ucds[name]
        else: ucd=None
        t.add_column(name,col,dtype=tipe,ucd=ucd)
        if name.lower()=='ra':
            print col[0]
            t.add_column('_Ra',[e.hours(s)*180/n.pi for s in col],dtype='<f8',ucd="pos.eq.ra;meta.main")    
        if name.lower()=='dec':
            t.add_column('_Dec',[e.degrees(s)*180/n.pi for s in col],dtype='<f8',ucd="pos.eq.dec;meta.main")    

#    t.add_column('Name',[s.split('\t')[0].strip() for s in lines if not s.startswith('#')],dtype=n.str)
#    t.add_column('Ra',[s.split('\t')[1].strip() for s in lines if not s.startswith('#')],dtype=n.str)
#    t.add_column('Dec',[s.split('\t')[2].strip() for s in lines if not s.startswith('#')],dtype=n.str)
#    t.add_column('S_nu_',[s.split('\t')[3].strip() for s in lines if not s.startswith('#')],dtype='<f8')
#    try: t.add_column('e_S_nu_',[s.split('\t')[4].strip() for s in lines if not s.startswith('#')],dtype='<f8')
#    except(IndexError): t.add_column('e_S_nu_',['0' for s in lines if not s.startswith('#')],dtype='<f8')
    if not opts.add is None:
        adds = opts.add.split(',')
        for add in adds:
            t.add_column(add.split(':')[0],[add.split(':')[1] for s in lines if not s.startswith('#')],dtype='S14')
    print file + ' > ' + '.'.join(file.split('.')[:-1])+'.vot'
    t.write('.'.join(file.split('.')[:-1])+'.vot')
    
