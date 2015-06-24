#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os
import glob

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--xx', default=None, help='Filenames for the xx pols.')
o.add_option('--yy', default=None, help='Filenames for the yy pols.')
o.add_option('--xy', default=None, help='Filenames for the xy pols.')
o.add_option('--yx', default=None, help='Filenames for the yx pols.')
opts,args = o.parse_args(sys.argv[1:])

if opts.xx == None and opts.yy == None and opts.xy == None and opts.yx == None:
    print 'This script is not for one file.'
    sys.exit()
elif opts.xx != None and opts.yy != None and opts.xy == None and opts.yx == None:
    xxfiles = glob.glob(opts.xx)
    yyfiles = glob.glob(opts.yy)
    files = zip(xxfiles,yyfiles)
    NPOLFILES = 2
elif opts.xx != None and opts.yy != None and opts.xy != None and opts.yx != None:
    xxfiles = glob.glob(opts.xx)
    yyfiles = glob.glob(opts.yy)
    xyfiles = glob.glob(opts.xy)
    yxfiles = glob.glob(opts.yx)
    files = zip(xxfiles,yyfiles,xyfiles,yxfiles)
    NPOLFILES = 4
else:
    print 'Invalid combination of linear pols... Exiting.'
    sys.exit()

print NPOLFILES

for filelist in files:
    outfile = filelist[0].replace('.xx.','.')
    print filelist, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dsum,dwgt = {}, {}
    curtime = None
    polcnt = 0
    for linpol in xrange(NPOLFILES):
        polcnt += 1
        uvi = a.miriad.UV(filelist[linpol])
        if polcnt == 1:
            uvo = a.pol.UV(outfile, status='new')
            uvo.init_from_uv(uvi)
        ants = a.scripting.parse_ants(opts.ant, uvi['nants'])
        a.scripting.uv_selector(uvi, ants)
        for p,d in uvi.all():
            pol=a.miriad.pol2str[uvi['pol']]
            uvo.write_pol(pol)
            uvo.write(p,d)
            
    uvo._wrhd('history',uvo['history'] + 'COMBINE_POL:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)


