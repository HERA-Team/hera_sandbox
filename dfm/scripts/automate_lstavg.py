#! /usr/bin/env python
import aipy as a
import numpy as np
import sys,os,ephem,optparse


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--lst_res',dest='dlst',type='float',default=30.,help='Resolution in seconds of output LST bins (default 30)')
o.add_option('--lst_rng', default='0_24', help='Range of LSTs to bin (in hours).')
o.add_option('--tfile',dest='tfile',type='float',default=600.,help='Length of time spanned by each input file.')
o.add_option('--tout',dest='tout',type='int',default=10,help='Length of output file (in sidereal minutes)')
opts,args = o.parse_args(sys.argv[1:])

RAstart,RAend = opts.lst_rng.split('_')
RAnow = RAstart
while( ephem.hours(RAnow) < ephem.hours(RAend) ):
    h,m = RAnow.split(':')
    h = int(h)
    m = int(m)+opts.tout
    if m >= 60: 
        h += 1
        m = m%60
    RAnext = '%d:%02d'%(h,m)
    this_lst_rng = RAnow+'_'+RAnext
    command = 'python add_lst_TS.py -C %s --lst_res=%f --lst_rng=%s --tfile=%d %s' \
        % (opts.cal,opts.dlst,this_lst_rng,opts.tfile,' '.join(args))
    print 'Working on LST range',this_lst_rng
    os.system(command)
    RAnow = RAnext
