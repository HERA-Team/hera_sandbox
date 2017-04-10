#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C, pylab as p
import optparse, sys, os
import glob

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--xx', default=None, help='Filenames for the xx pols.')
o.add_option('--yy', default=None, help='Filenames for the yy pols.')
o.add_option('--xy', default=None, help='Filenames for the xy pols.')
o.add_option('--yx', default=None, help='Filenames for the yx pols.')
opts,args = o.parse_args(sys.argv[1:])

xxfiles = glob.glob(opts.xx)
yyfiles = glob.glob(opts.yy)
xyfiles = glob.glob(opts.xy)
yxfiles = glob.glob(opts.yx)
files = zip(xxfiles,yyfiles,xyfiles,yxfiles)
linpols = ['xx','yy','xy','yx']

for filelist in files:
    outfile = os.path.commonprefix(filelist) + os.path.basename(filelist[0]).split('.')[-1]
    print filelist, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dlist = {}
    curtime = None
    for polcnt, linpol in enumerate(linpols):
        dlist[linpol] = {}
        uvi = a.miriad.UV(filelist[polcnt])
        ants = a.scripting.parse_ants(opts.ant, uvi['nants'])
        a.scripting.uv_selector(uvi, ants)
        #uvi.select('decimate', opts.decimate, opts.decphs)
        for (crd,t,(i,j)),d in uvi.all():
            bl = a.miriad.ij2bl(i,j)
            if t != curtime:
                curtime = t
                dlist[linpol][t] = {}
            dlist[linpol][t][bl] = d
        uvi.rewind()

    print '    Writing output file'
    uvo = a.pol.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'pol':a.miriad.str2pol['xx']})
    for pol in dlist.keys():
        uvo.write_pol(pol)
        times = n.array(dlist[pol].keys()).astype(float)
        times.sort()
        for t in times:
            for bl in dlist[pol][t].keys():
                i,j = a.miriad.bl2ij(bl)
                p = (n.array([0.,0.,0.]),t,(i,j))
                d = dlist[pol][t][bl]
                uvo.write(p,d)
            
    uvo._wrhd('history',uvo['history'] + 'MERGE_POL:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




