#! /usr/bin/env python
'''
Takes outputs from 21cmsense and finds P(k) values for the noise
'''

def r(number):
    return '%5.0f'%n.round(number,1)
    #return number

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['xtick.labelsize']= 18
mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['lines.markersize']= 12

mpl.rcParams['legend.numpoints']  = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['figure.dpi'] = 500
mpl.rcParams['savefig.dpi'] = 500
mpl.rcParams['savefig.format'] ='png'

from pylab import *
import numpy as n,sys,os,re
from capo.cosmo_units import *
from capo import pspec 
import optparse,glob
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from capo import pspec,eor_results,cosmo_units
import ipdb
import capo
import matplotlib.pyplot as p
from IPython import embed

o = optparse.OptionParser()
o.set_usage('21cmsense_to_pk.py [options]')
o.set_description(__doc__)
o.add_option('--time',default=0,type=float,
    help='integration time in hours')
o.add_option('--models', type='str',
        help='a string that globs a list of 21cmfast pspec output files')
o.add_option('--freq', type='float',default=150,
        help='Supply frequency for output P(k) value (default=150MHz)')
opts,args = o.parse_args(sys.argv[1:])

noisefiles = sort(args)
print "loading a few noise model redshifts and building a 2d (z,k) cubic fit model"
#re_f = re.compile(r'paper_dcj_lstcnt_sa_pess_(\d+\.\d+)\.npz')
re_f = re.compile('psa6240_v003drift_mod_(\d+\.\d+)\.npz')#glob should load them in increasing freq order

noises = []
noise_ks = []
noise_freqs = []
nk_grid = n.linspace(0,1)*0.5+0.01
for noisefile in noisefiles:
    noise = n.load(noisefile)['T_errs']
    noise_k = n.load(noisefile)['ks']

    #look for infs at low k and re-assign large values
    #noise[ n.logical_and(n.isinf(noise), noise_k < .15)] = 1e+6

    bad = n.logical_or(n.isinf(noise),n.isnan(noise))
    noise = noise[n.logical_not(bad)] 
    noise_k = noise_k[n.logical_not(bad)]
    #keep only the points that fall in our desired k range
    good_k = noise_k < nk_grid.max()
    noise = noise[noise_k<nk_grid.max()]
    noise_k = noise_k[noise_k<nk_grid.max()]
    print noisefile,n.max(noise),

    noise_k = n.insert(noise_k,0,0)
    noise = n.insert(noise,0, 0)

    #nk3 = noise_k**3/(2*np.pi**2)
    #tmp = n.linalg.lstsq( n.vstack([nk3,n.zeros((3,len(nk3)))]).T,noise)
    tmp = n.polyfit(noise_k,noise,3,full=True)
    noise = n.poly1d(tmp[0])(nk_grid)
    noises.append(noise)
    noise_ks.append(nk_grid)
    f = float(re_f.match(noisefile).groups()[0])*1e3 #sensitivity freq in MHz
    print f
    noise_freqs.append(f) 
    #embed()

noise_k_range = [n.min(n.concatenate(noise_ks)),n.max(n.concatenate(noise_ks))]
nk_count = n.mean([len(myks) for myks in noise_ks])
nks = n.linspace(noise_k_range[0],noise_k_range[1],num=nk_count)
noise_interp = n.array([interp(nks,noise_ks[i],noises[i]) for i in range(len(noises))])
NK,NF = n.meshgrid(nks,noise_freqs)
#noise_freqs = n.array(noise_freqs)
POBER_NOISE = interp2d(NK,NF,noise_interp,kind='linear')#build 2d interpolation model


noise_delta2= POBER_NOISE(nk_grid, opts.freq)

k3=nk_grid**3/(2*np.pi**2)
#pk = n.linalg.lstsq( n.vstack([k3,n.zeros((3,len(k3)))]).T,noise_delta2)[0][0]
pk = noise_delta2/k3

print 'Noise P(k) [mk^2]/[hMpc^-1]^3'
print '\tmean: {0:.4e}'.format( n.mean(pk))
print '\tmedian: {0:.4e}'.format( n.median(pk))
