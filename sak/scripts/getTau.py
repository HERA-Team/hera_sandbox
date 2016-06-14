#! /usr/bin/env python
"""
Calculates the value of tau such that the statistic:

X^2(nu,t) = sum_{baselines} V^xy_b - V^{yx}_b exp(-2 pi i nu tau)

is minimized across the array. 
Assumes data naming convention zen.2451234.56789.pp.uv 
and that the directory that holds xy data also holds yx

Returns npz of the chisquared array and optionally applies the
best chisquare to the V^{yx} data

Author: Saul Aryeh Kohn
"""

import aipy, os, sys, glob, numpy as np, optparse, capo
from matplotlib import pylab

o = optparse.OptionParser()
o.set_usage('getTau.py [options] *xy.uv')
o.set_description(__doc__)
o.add_option('--taumax',dest='taumax',default=30.,type='float',help='delays -taumax<=t<=taumax will be investigated. Units: ns')
o.add_option('--taures',dest='taures',default=0.01,type='float',help='resolution of delay sweep. Units: ns')
o.add_option('-A','--apply',dest='app',default=False,action='store_true',help='Create new yx file with optimal tau applied')
opts,args = o.parse_args(sys.argv[1:])

nu_arr = np.linspace(.1, .2, num=203) #GHz
tau_arr = np.arange(-1*opts.taumax, opts.taumax, opts.taures) #nanseconds
angle_arr = np.exp(-2*np.pi*1.j*np.outer(nu_arr,tau_arr))

for ii,f in enumerate(args):
    xyfile = f
    yxfile = np.array(f.split('.'))
    yxfile[3] = 'yx'
    yxfile = '.'.join(yxfile)
    
    _,dxy,fxy = capo.arp.get_dict_of_uv_data([xyfile], antstr='cross', polstr='xy') #dxy[(6,32)]['xy'].shape = (14,203)
    _,dyx,fyx = capo.arp.get_dict_of_uv_data([yxfile], antstr='cross', polstr='yx')
    
    chisquared = np.zeros( (tau_arr.shape[0], dxy[dxy.keys()[0]]['xy'].shape[0], nu_arr.shape[0]) ) #tau x time x freq
    
    for jj in range(angle_arr.shape[1]):
        calc = 0.
        for bl in dxy.keys():
            xy = dxy[bl]['xy']
            yx = dyx[bl]['yx']
            calc += np.power(np.absolute(xy - yx*angle_arr[:,jj]),2.)
        chisquared[jj,:,:]=calc
    
    tf = yxfile.split('.')
    tf[0] = 'tau'
    tf = '.'.join(tf)+'.npz'
    print 'Saving chisq %s'%tf
    np.savez(tf, chisq=chisquared, tlim=opts.taumax, tres=opts.taures)
    del(dxy);del(dyx);del(fxy);del(fyx)
    
    if opts.app:
        uvofile = yxfile+'T'
        print yxfile,'->',uvofile
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
            
        uvi = aipy.miriad.UV(yxfile)
        uvo = aipy.miriad.UV(uvofile, status='new')
        
        mins = np.argmin(chisquared,axis=0)
        final_taus = np.zeros(mins.shape)
        for t in range(mins.shape[0]): final_taus[t,:] = tau_arr[mins[t,:]] #select-out taus per freq
        x_angle = np.exp(-2*np.pi*1.j*nu_arr*final_taus)
        angles = {'xA':x_angle} #have to 'hide' this in a dict to be able to pass it into the mfunc
        
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
        uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist="APPLY_TAUCALC: yx phased to minimize V according to %s \n"%tf)
        del(uvo); del(curtime); del(t_index); del(angles)
    else:
        continue



