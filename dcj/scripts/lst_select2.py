#!/usr/global/paper/bin/python
"""
Read in a list of files and outputs file list with LSTs within RA range.
"""

import aipy as a, sys, optparse, ephem,string,numpy as n

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
o.add_option('--lst_pad',dest='lst_pad',default=0,
    help="LST search pad.  Increase the search range by this many hours.")
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .1, .1, 1)
sun = ephem.Sun()
if len(args) == 0: args = [a.phs.ephem2juldate(ephem.now())]
if not opts.src is None: 
    srclist,cutoff = a.scripting.parse_srcs(opts.src)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff)
    cat.compute(aa)
else: srclist =[]
#parse the ra range
if not opts.ra_rng is None:
    ra1,ra2 = map(ephem.hours, opts.ra_rng.split('_'))
    if ra1>ra2: ra1,ra2 = ra2,ra2
else:
    ra1,ra2  = (0,0)
active=0
pad=ephem.hours(opts.lst_pad)
for s in args:
    jd = float(string.join(s.split('.')[1:3],'.'))
    aa.set_jultime(jd)
    #print s,aa.sidereal_time(),repr(aa.sidereal_time())
    sun.compute(aa)
    #print ra1,aa.sidereal_time(),ra2
    if aa.sidereal_time()>(ra1-pad) and aa.sidereal_time()<(ra2+pad):
        if active==0 and not opts.dchar is None: print opts.dchar
        print s
        active=1
        if not opts.debug is None: print 'ra range'
    elif aa.sidereal_time()>(sun.ra-20.0*n.pi/180-pad) and aa.sidereal_time()<(sun.ra+20*n.pi/180+pad) and opts.sun:
        if active==0 and not opts.dchar is None: print opts.dchar
        print s
        active =1
        if not opts.debug is None: print 'sun'
    elif len(srclist)>0:
        for src,obj in cat.iteritems():
            if aa.sidereal_time()>(obj.ra-20.0*n.pi/180-pad) and aa.sidereal_time()<(obj.ra+20*n.pi/180+pad):
                if active==0 and not opts.dchar is None: print opts.dchar
                print s
                active =1
                if not opts.debug is None: print 'src'
                break
    else:active=0
    if not opts.debug is None: print aa.sidereal_time()
#    sun.compute(aa)
#    print 'LST:', aa.sidereal_time(),
#    print '     Julian Date:', jd,
#    print '     Day:', aa.date
#    print 'Sun is at (RA, DEC):', str((str(sun.ra), str(sun.dec)))
#    print '-' * 70

