#! /usr/bin/env python
import numpy as n, pylab as p, sys
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse, sys

o = optparse.OptionParser()
o.add_option('--output', type='string', default=None,
    help='Output filename to save power spectrum values.')
opts,args = o.parse_args(sys.argv[1:])


### GETTING SIGLOSS DATA ###

pCs,pIs,pCvs,pIvs = [],[],[],[]
pCs_err,pIs_err = [],[]
freq = []
for inject in glob.glob('inject_*'):
    print 'Reading', inject
    pspecs = glob.glob(inject + '/pspec_boot*.npz')
    inject = float(inject.split('_')[-1])
    pC_avg, pI_avg, pCv_avg, pIv_avg = [], [], [], []
    pC_spec, pI_spec = [], []
    pcf,pif = [],[]
    for pspec in pspecs:
        npz = n.load(pspec)
        try: freq = npz['freq']
        except: pass
        pC,pI = npz['pk_vs_t'], npz['nocov_vs_t']
        try: pCv,pIv =  npz['pCv'], npz['pIv'] #(#chan, #times)
        except: pass
        kpls = n.array(npz['kpl'])
        pC_avg.append(n.average(pC.real)) #avg over freq and time
        pcf.append(pC.real)
        pI_avg.append(n.average(pI.real))
        try: 
            pCv_avg.append(n.average(pCv.real,axis=1)) #spectrum
            pIv_avg.append(n.average(pIv.real,axis=1))
        except: pass
        #pC_spec.append(n.average(pC.real, axis=1))
        #pI_spec.append(n.average(pI.real, axis=1))
        #pC_spec.append(pC.real)
        #pI_spec.append(pI.real)
    pCs.append(n.average(pC_avg)) 
    pCs_err.append(n.std(pC_avg)/n.sqrt(len(pspecs)))
    pIs.append(n.average(pI_avg))
    pIs_err.append(n.std(pI_avg)/n.sqrt(len(pspecs)))
    try:  
        pCvs.append(n.average(pCv_avg,axis=0)) #spectrum
        pIvs.append(n.average(pIv_avg,axis=0))
    except: pass
    #pC_spec = n.average(pC_spec, axis=0)
    #print inject, pspec, pC_avg, pI_avg
    #p.figure(1)
    #p.subplot(211); p.loglog([pI_avg], [pC_avg], 'k.')
    #p.subplot(212); p.loglog([pC_avg], [pI_avg/pC_avg], 'k.')

pIs,pCs,pCvs,pIvs = n.array(pIs), n.array(pCs), n.array(n.average(pCvs,axis=0)), n.array(n.average(pIvs,axis=0)) #avg over inject #s
pIs_err,pCs_err = n.array(pIs_err), n.array(pCs_err)

#Plot pCv points from pspec_cov_v004 output, not sigloss output
pCv_points = n.load('pspec_C.npz')['pk']
pCvs = n.abs(pCv_points) #XXX overwrites pCvs from inject outputs
pIv_points = n.load('pspec_I.npz')['pk']
pIvs = n.abs(pIv_points)

###Build an interpolator to find sigloss factors###
sig_factor_interp = interp1d(n.abs(pCs_full.ravel()), n.abs(pIs_full.ravel())/n.abs(pCs_full.ravel()),
                        kind='linear',bounds_error=False,fill_value=0)
order = n.argsort(n.abs(pCs)) #set up interpolator to work even if pCs are out of order
pCs_order = n.abs(pCs[order])
pIs_order = n.abs(pIs[order])
try: sig_factor_interp = interp1d(pCs_order, pIs_order/pCs_order,kind='linear',bounds_error=False,fill_value=0)
except: pass

p.plot(n.abs(kpls),n.abs(pCvs),'k.')
#p.errorbar(kpls, n.abs(pCvs), yerr=2*pCvs_err, fmt='k.', capsize=0)
p.grid()

try:
    for k,kpl in enumerate(kpls):
        print ("%5.5f" % kpl), ':', ("%5.5f" % n.abs(pCvs[k])), ':', float(sig_factor_interp(n.abs(pCvs[k])))
    maxfactor = float(sig_factor_interp(n.max(n.abs(pCvs))))
    print "Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,maxfactor) #n.max(sig_factors))
except: pass

p.show()

if opts.output != None:
    print 'Writing '+opts.output+'.npz'
    n.savez(opts.output+'.npz', kpl=kpls, pCv=pCvs, pIv=pIvs, factor=maxfactor, cmd=' '.join(sys.argv))

#p.show()
p.close()
