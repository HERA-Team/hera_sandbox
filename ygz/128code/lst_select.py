#!/usr/bin/env python
"""
Read in a list of files and outputs file list with LSTs within RA range.
"""

import aipy as a, sys, optparse, ephem,string,numpy as n,logging,os

o = optparse.OptionParser()
o.set_usage('lst_select [options] <files>')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('--ra',dest='ra_rng',
    help='A range RA1_RA2 of right-ascensions to select for.  Default: None.')
o.add_option('--sun',dest='sun',action='store_true',
    help="Find files near local noon")
o.add_option('-d',dest='dchar',type='str',
    help="""Divide list into seperate days with this character.""")
o.add_option('--debug',dest='debug',action='store_true',
    help="Output optional debug info.")
o.add_option('--lst_list',dest='lst_list',action='store_true',
    help="Output optional lst info.")
o.add_option('--lst_pad',dest='lst_pad',default=0,type='float',
    help="LST search pad.  Increase the search range by this many hours.")
o.add_option('--suntime',default='e',
    help='Sun up = y, sun down = n, either = e [default=e]')
opts, args = o.parse_args(sys.argv[1:])

if opts.debug:
    logging.basicConfig(level=logging.DEBUG)
else: 
    logging.disable(logging.WARNING)
log = logging.getLogger('lst_select')

aa = a.cal.get_aa(opts.cal, .1, .1, 1)
sun = ephem.Sun()
if len(args) == 0: args = [a.phs.ephem2juldate(ephem.now())]
if not opts.src is None:
    srclist,cutoff,cats = a.scripting.parse_srcs(opts.src,opts.cat)
    if not opts.cal is None:
        cat = a.cal.get_catalog(opts.cal,srcs=srclist, cutoff=cutoff,catalogs=cats)
    else:
        cat = a.src.get_catalog(srcs=srclist, cutoff=cutoff,catalogs=cats)
    cat.compute(aa)
else: srclist =[]
#parse the ra range
if not opts.ra_rng is None:
    ra1,ra2 = map(ephem.hours, opts.ra_rng.split('_'))
    if ra1>ra2: ra1,ra2 = ra2,ra1
else:
    ra1,ra2  = (0,0)
active=0
is_listed = False
pad=ephem.hours(opts.lst_pad*a.img.deg2rad)
for s in args:
    #print os.path.basename(s).split('.')
    jd = float(string.join(os.path.basename(s).split('.')[1:3],'.'))
    aa.set_jultime(jd)
    #print s,aa.sidereal_time(),repr(aa.sidereal_time())
    sun.compute(aa)
    #print ra1,aa.sidereal_time(),ra2
    if aa.sidereal_time()>(ra1-pad) and aa.sidereal_time()<(ra2+pad):
        if active==0 and not opts.dchar is None: print opts.dchar
#        print s,
        active=1
        is_listed=True
        if not opts.debug is None: print 'ra range',
    elif aa.sidereal_time()>(sun.ra-20.0*n.pi/180-pad) and aa.sidereal_time()<(sun.ra+20*n.pi/180+pad) and opts.sun:
        if active==0 and not opts.dchar is None: print opts.dchar
#        print s,
        active =1
        is_listed=True
        if not opts.debug is None: print 'sun',
    elif len(srclist)>0:
        for src,obj in cat.iteritems():
            if opts.debug: print aa.sidereal_time(),obj.ra
            if aa.sidereal_time()>(obj.ra-2.*n.pi/12-pad) and aa.sidereal_time()<(obj.ra+1.*n.pi/12+pad):
                if active==0 and not opts.dchar is None: print opts.dchar
#                print s,
                active =1
                is_listed = True
                if not opts.debug is None: print obj.src_name,
                break
    else:
        active=0
        is_listed = False
    if opts.suntime=='n' and sun.alt>0:
        active=0
        is_listed=False
    elif opts.suntime=='y' and sun.alt<0:
        active=0
        is_listed=False
    if is_listed: print s
    if not opts.lst_list is None and is_listed: print "\t",aa.sidereal_time()
    if not opts.debug is None and opts.lst_list is None and is_listed: print "."
    is_listed=False
#    elif not opts.lst_list is None and not is_listed: print
#    sun.compute(aa)
#    print 'LST:', aa.sidereal_time(),
#    print '     Julian Date:', jd,
#    print '     Day:', aa.date
#    print 'Sun is at (RA, DEC):', str((str(sun.ra), str(sun.dec)))
#    print '-' * 70

