#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import atpy

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

rfi = n.load('/data3/paper/eor_reduce/FlagSummary.npz')

aa = a.cal.get_aa(opts.cal,n.array([.150]))

lst = []
ofilename = 'rfi_flags.vot'
if True:
    t = atpy.Table()
    t.add_column('Aflgs', rfi['allflags'], dtype='<f8')
    t.add_column('Mflgs', rfi['manflags'], dtype='<f8')
    t.add_column('Cflgs', rfi['corrflags'], dtype='<f8')
    t.add_column('Sflgs', rfi['statflags'], dtype='<f8')
    jds = rfi['JD']
    for jd in jds:
        aa.set_jultime(jd)
        lst.append(aa.sidereal_time()*12/n.pi)
    t.add_column('jds', jds, dtype='<f8')
    t.add_column('lsts', lst, dtype='<f8')
    t.write(ofilename,overwrite=True)
