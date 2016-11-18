#! /usr/bin/env python

"""
This calculates the average of XY/YX ratios over all baselines in the uv file
in order to estimate the x-to-y phase offset (minimizing Stokes V):

exp(-2pi i nu <tau_xy>) = 1/N_bl \sum^{N_bl}_{i,j} Vxy_{ij}/Vyx_{i,j}

"""


import aipy, os, sys, glob, numpy as np, optparse, capo
from matplotlib import pylab

o = optparse.OptionParser()
o.set_usage('getTau_v2.py [options] *xy.uv')
o.set_description(__doc__)
#o.add_option('--taumax',dest='taumax',default=30.,type='float',help='delays -taumax<=t<=taumax will be investigated. Units: ns')
#o.add_option('--taures',dest='taures',default=0.01,type='float',help='resolution of delay sweep. Units: ns')
o.add_option('-A','--apply',dest='app',default=False,action='store_true',help='Create new yx file with optimal tau applied')
opts,args = o.parse_args(sys.argv[1:])

nu_arr = np.linspace(.1, .2, num=203) #GHz

for ii,f in enumerate(args):
    xyfile = f
    yxfile = np.array(f.split('.'))
    yxfile[3] = 'yx'
    yxfile = '.'.join(yxfile)
    
    print 'Reading %s, %s'%(xyfile, yxfile),
    tin,dxy,fxy = capo.arp.get_dict_of_uv_data([xyfile], antstr='cross', polstr='xy') #dxy[(6,32)]['xy'].shape = (14,203)
    _,dyx,fyx = capo.arp.get_dict_of_uv_data([yxfile], antstr='cross', polstr='yx')
    print 'Done'
    
    """
    I'm going to calculate the ratios, and then sum those ratios across the BLs.
    Is this the better order for reducing noise?
    XXX Do error propagation problem
    """
    rdict = {}
    for bl in dxy.keys(): rdict[bl] = dxy[bl]['xy']/dyx[bl]['yx'] #ratio
    sumarr = np.zeros_like(dxy[dxy.keys()[0]]['xy'])
    for bl in rdict.keys(): sumarr+=rdict[bl] #sum
    sumarr = sumarr/len(rdict.keys()) #normalize
    
    tf = yxfile.split('.')
    tf[0] = 'tau'
    tf = '.'.join(tf)+'.npz'
    print '    Saving %s'%tf
    np.savez(tf,angle=sumarr)
    
    
    if opts.app:
        uvofile = yxfile+'T'
        print yxfile,'->',uvofile
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
            
        uvi = aipy.miriad.UV(yxfile)
        uvo = aipy.miriad.UV(uvofile, status='new')
        angles = {'xA':sumarr} #have to 'hide' this in a dict to be able to pass it into the mfunc
        curtime=0
        t_index=-1
        def mfunc(uv,p,d,f):
            global curtime
            global t_index
            global angles
            uvw, t, (i,j) = p
            if t!=curtime:
                t_index+=1
                curtime =t
            angle = angles['xA'][t_index,:]
            d = d*angle
            return p,d,f
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist="getTau_v2: yx phased to minimize V w.r.t. %s"%xyfile)
        del(uvo); del(curtime); del(t_index); del(angles)
        del(uvi)
        uvi = aipy.miriad.UV(xyfile)
        uvo = aipy.miriad.UV(xyfile+'T', status='new') #makes next processing stages easier if they have the same ending.
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,append2hist="getTau_v2: yx phased to minimize V w.r.t this file")
        del(uvo)
        del(uvi)
    else:
        continue
