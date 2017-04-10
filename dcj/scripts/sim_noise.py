#! /usr/bin/env python
"""
Take a dataset and replace it with noise
"""
import numpy as np, aipy as a, optparse, os, sys, ephem,capo

o = optparse.OptionParser()
o.set_usage('sim_noise.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise (Jy) added to each UV sample of simulation.')
o.add_option('--Trcvr',type=float,
    help="""Flat spectrum noise level in Kelvins. 
    When set triggers calculation of sky-model+Trcvr noise and overrides
    --noiselev.""")
o.add_option('--usevar',action='store_true',
    help='use variance output by lstbinner for noise spectrum. overrides -n but not Trcvr')
o.add_option('--inttime',type=float,
    help='override the integration time')
o.add_option('--outdir',type=str,
    help='write files to this directory')
opts, args = o.parse_args(sys.argv[1:])


def complex_noise(size,noiselev):
    if noiselev<=0 or np.isinf(noiselev):return np.zeros(size)
    noise_real = np.random.normal(size=size, scale = noiselev)/np.sqrt(2)
    noise_imag = np.random.normal(size=size, scale = noiselev)/np.sqrt(2)
    noise = noise_real +1j*noise_imag
    return noise

curtime=None
if not opts.Trcvr is None:
    def noise(uv,p,d,f):
        #here we assume a constant sky temperature of 180K at 180MHz with a spectral index of -2.88. 
        # adding in LST dependence is a good task for the interested student
        # lst = uv['lst']
        df = uv['sdf'] * 1e9
        ndays = uv['cnt']
        Tsys = 180*(freqs/0.18)**-2.55 + opts.Trcvr #system temp in K
        Tsys *= 1e3 #system temp in mK
        if not opts.inttime is None: inttime=opts.inttime
        else:inttime=uv['inttime']
        Trms = Tsys/np.sqrt(df*inttime * ndays)
        Vrms = Trms/jy2T  #jy2T is in units of mK/Jy
        noise = np.array([complex_noise(1,v) for v in Vrms])
        noise.shape = d.shape
        return p, noise, f 
elif opts.usevar:#EXPERIMENTAL
    def noise(uv,p,d,f):
        Vrms = np.sqrt(uv['var'].real)
        noise = np.array([complex_noise(1,v) for v in Vrms])
        noise.shape = d.shape
        return p,noise,f
else:
    def noise(uv, p, d, f):
        # Add gaussian noise.
        d = complex_noise(d.size,opts.noiselev)# * aa.passband(i,j)
        return p, np.where(f, 0, d), f    

for filename in args:
    if not opts.outdir is None:
        uvofile = opts.outdir+'/'+os.path.basename(filename)+'s'
    else:
        uvofile = filename + 's'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    freqs = a.cal.get_freqs(uvi['sdf'], uvi['sfreq'], uvi['nchan'])
    jy2T = capo.pspec.jy2T(freqs)
    #a.scripting.uv_selector(uvi, opts.ant)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=noise, raw=True, append2hist="sim_noise: " + ' '.join(sys.argv) + '\n')
