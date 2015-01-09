#! /usr/bin/env python

import aipy as a, numpy as n
import optparse, sys, os
from scipy.io import readsav
import pylab as plt

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True,src=True)
o.add_option('--mirdir',default='',help='Directory where corresponding miriad files live.')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(opts.mirdir+os.path.basename(args[0]).split('_')[0])
aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
del(uv)

for filename in args:

    #read the sav file
    savfile = readsav(filename)
    pol = filename.split('_')[-1][:-4].lower()
    try: vis_array = savfile['vis_ptr']
    except(KeyError): vis_array = savfile['vis_model_ptr']
    obs = savfile['obs']
    times = obs['baseline_info'][0]['JDATE'][0]
    ants1 = obs['baseline_info'][0]['TILE_A'][0]
    ants2 = obs['baseline_info'][0]['TILE_B'][0]
    ntimes = len(times)
    nbls = obs['NBASELINES'][0]

    #reorder vis_array
    time_order = n.argsort(times)
    times = times[time_order]
    _vis_array = n.zeros_like(vis_array)
    for ntime in xrange(ntimes):
        correct_time = time_order[ntime]
        _vis_array[ntime*nbls:(ntime+1)*nbls] = vis_array[correct_time*nbls:(correct_time+1)*nbls]
    vis_array = _vis_array

    #create matching uv
    fhdname = os.path.basename(filename)
    mirname = fhdname.split('_')[0]
    if opts.mirdir: mirpath = opts.mirdir + mirname
    else: mirpath = mirname
    print 'reading UV file', mirname
    if os.path.exists(mirpath+'H'+pol):
        print '\tFile exists... skipping.'
        continue
    uv = a.pol.UV(mirpath)
    uvo = a.pol.UV(mirpath+'H'+pol, status='new')
    uvo.init_from_uv(uv)

    #read time for unphasing
    uvi = a.miriad.UV(mirpath[:-1])
    (uvw,t0,(i,j)), dat = uvi.read()
    del(uvi)
    #t0 = float(mirname.split('.')[1]+'.'+mirname.split('.')[2])
    print 'initial time', t0
    aa.set_jultime(t0+5.*60./a.const.s_per_day)
    RA = str(aa.sidereal_time())
    dec= str(aa.lat)
    opts.src = RA+'_'+dec
    del(uv)
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]

    #loop through and write
    timecnt,blcnt = 0,0
    aa.set_active_pol(pol)
    for ind, d in enumerate(vis_array):
        #calculate the right timestamp
        if blcnt == nbls:
            blcnt = 0
            timecnt += 1
        else:
            blcnt += 1
        time = times[timecnt]
        
        i = ants1[ind] - 1 
        j = ants2[ind] - 1       

        p = (n.array([0.,0.,0.]),time,(i,j))

        #unphase the data
        aa.set_jultime(time)
        src.compute(aa)
        print src
        d = n.conj(d) #you need to conjugate all the data and i dont know why
        d = aa.unphs2src(d,src,i,j)
 
        f = n.zeros_like(d) + n.where(d == 0.,1,0)
        uvo.write_pol(pol)
        uvo.write(p,d,f)
    del(src)

    uvo._wrhd('history', uvo['history'] + 'sav2miriad: created file from %s' % fhdname + '\n')
    del(uvo)
