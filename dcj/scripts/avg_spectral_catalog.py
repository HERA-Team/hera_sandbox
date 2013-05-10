#! /usr/bin/env python
"""
Average the output of beam_src_vs_ha into a single txt file.
"""
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n,atpy
import optparse, sys, scipy.optimize
import capo as C
from pylab import *
from scipy import optimize

CAT=True
#output='psa64_pic_stripe.txt'
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
o.add_option('--freqs',type='str',default='100_200_10',
    help='spectral range in MHz start_stop_chan default=100_200_10')
opts,args = o.parse_args(sys.argv[1:])

def filename2src(f):
    #return f.split('_')[-1]
    return f.split('_')[0].strip()
def average_spectrum(x,y,bins):
    myspectrum = []
    myerrors = []
    myfreqs = []
    for i in range(1,len(bins)):
        myspectrum.append(n.mean(y[n.logical_and(x>bins[i-1],x<bins[i])]))
        myerrors.append(n.std(y[n.logical_and(x>bins[i-1],x<bins[i])]))
        myfreqs.append((bins[i]+bins[i-1])/2)
    myspectrum=n.array(myspectrum)
    myfreqs = n.array(myfreqs)
    myerrors = n.array(myerrors)
    return myfreqs,myspectrum,myerrors


freqbins = map(float,opts.freqs.split('_'))
freqbins = n.arange(freqbins[0],freqbins[1],freqbins[2])
for i,filename in enumerate(args):
    srcname=filename2src(filename)
    D = n.load(filename)
    freq = D['freq']*1e3
    spec = n.ma.masked_invalid(D['spec'])
    psa64freqs,psa64fluxes,psa64errors = average_spectrum(freq,n.real(spec),freqbins)
    if i==0: print "#FREQS[MHz]=%s"%(','.join(map(str,psa64freqs)))
    print "%s\t%s\t%s"%(srcname,
        ','.join(map(str,n.round(n.real(psa64fluxes),2))),
        ','.join(map(str,n.round(psa64errors,3))))


