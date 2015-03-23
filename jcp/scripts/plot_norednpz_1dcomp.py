#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os
from matplotlib import rc
rc('text',usetex=True)
rc('font', size=22)

o = optparse.OptionParser()
o.add_option('--nobin', action='store_true',
    help='Do not bin by u-magnitude')
o.add_option('--umax', type='float', default=n.Inf,
    help='Only show baselines shorter than this value')   
o.add_option('--umin', type='float', default=0.,
    help='Only show baselines longer than this value')   
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

fq = .1525
aa = a.cal.get_aa(opts.cal, n.array([fq]))

uDat, uWgt = {}, {}
for npzfile in args:
    print 'Reading...', npzfile
    dat = n.load(npzfile)
    kpl = dat['kpl']
    keys = dat.files[:]
    keys.remove('kpl')
    
    #bin in umag
    for bl in keys:
        if bl[0] == 'w': continue
        wbl = 'w'+str(bl)
        i,j = a.miriad.bl2ij(bl)
        #if i != 49 and j != 49: continue
        if i == 40 or j == 40: continue
        if i == 55 or j == 55: continue
        #if i == 59 or j == 59: continue
        crd = aa.get_baseline(i,j)*fq
        umag = (crd[0]**2 + crd[1]**2)**.5
        #pick predominantly east-west baselines
        #if crd[0]**2 < .5 * umag**2: continue
        if umag > opts.umax: continue
        if umag < opts.umin: continue
        if opts.nobin: umag = bl
        else:
            umag = str(2**int(n.ceil(n.log2(umag.clip(0.5,n.Inf)))))
            #umag = str(2**int(n.around(n.log2(umag.clip(0.5,n.Inf)))))
        uDat[umag] = uDat.get(umag,0) + dat[bl]
        uWgt[umag] = uWgt.get(umag,0) + dat[wbl]
    
keys = uDat.keys()
for umag in keys:
    uDat[umag] /= uWgt[umag]

keys = n.array(keys).astype(float)    
keys.sort()

colors = ['k','k','blue','green','red','cyan','purple','k']
labels = ['k','k','4-8','8-16','16-32','32-64','64-12','k']
for color,umag in enumerate(keys):
    umag = str(umag).split('.')[0]
    if umag in ['1','2','4','256']: continue #CLUDGE!
    f = n.abs(kpl)**3*(2*n.pi**2)**-1
    if opts.nobin: label = str(a.miriad.bl2ij(umag))
    #else: label = str(umag)
    else: label = labels[color]
    #log-bin in k
    if True:
        #binstart = 0.454
        binstart = [0,0,.123,.144,.185,.283,.454,0][color]
        nbins = 9
        logks = n.logspace(n.log10(binstart),1,nbins)
        binf = logks**3 * (2*n.pi**2)**-1
        binned_data, bin_counts = n.zeros(nbins),n.zeros(nbins)
        inds1 = n.digitize(n.abs(kpl),logks)
        inds2 = n.where(n.abs(kpl > binstart))[0]
        inds = inds1[inds2]
        #print inds1, inds2, inds
        for ind,dat in zip(inds,n.real(uDat[umag])[inds2]):
            binned_data[ind] += dat
            bin_counts[ind] += 1
        binned_data[bin_counts > 0] /= bin_counts[bin_counts > 0]
        #select unbinned data
        inds = n.where(n.abs(kpl) < binstart)
        ks,pspec = kpl.take(inds)[0],n.real(uDat[umag]).take(inds)[0]
        #average positive and negative ks
        pks = n.where(ks > 0)[0]
        nks = n.where(ks < 0)[0]
        new_pspec = n.zeros_like(pks)
        for ind,pind in enumerate(pks):
            nind = n.where(-1*ks == ks[pind])[0]
            #print pind,ks[pind],pspec[pind],nind,ks[nind],pspec[nind]
            new_pspec[ind] = (pspec[pind] + pspec[nind])/2
        ks,pspec = ks[pks],new_pspec
        #print umag, ks,pspec
        f = n.abs(ks)**3 * (2*n.pi**2)**-1
        p.loglog(n.abs(ks),2*f*n.abs(pspec),lw=2,marker='o',label=label,color=colors[color])
        p.loglog(logks,(2*binf*n.abs(binned_data)),lw=2,marker='o',color=colors[color])
        p.loglog([n.abs(ks)[-1],logks[1]],[(2*f*n.abs(pspec))[-1],(2*binf*n.abs(binned_data))[1]],color=colors[color],lw=2)
        print n.abs(ks)[-1], (f*n.abs(pspec))[-1]
        print logks[1], (binf*n.abs(binned_data))[1]
    #p.loglog(logks,f*n.abs(binned_data),label=label)
    #p.loglog(n.abs(kpl),f*n.abs(n.real(uDat[umag])),label=label)
    #p.loglog(n.abs(kpl),f*n.abs(n.imag(uDat[umag])),label=label)

    #p.loglog([binstart,binstart],[1,1e18],color=colors[color])
filename = '/data3/paper/arp/delay_pspec_paper/plots/fg_vs_umag_vs_fq_nonoise.npz'
npz = n.load(filename)
# Overlay foreground spectra from simulations
colors = ['green','red','cyan','purple']
inds = [0,1,2,3]
for ind, ks, fqs, specs in zip(inds, npz['ks'], npz['fqs'], npz['spec']):
    valid = n.where(n.abs(fqs - fq) < .002, 1, 0)
    valid *= n.where(specs == 1e6, 0, 1)
    #print valid.shape, ks.shape, fqs.shape, specs.shape, valid.sum()
    valid = valid.flatten()
    ks = n.abs(ks.flatten().compress(valid))
    specs = n.abs(specs.flatten().compress(valid))
    #print ks
    #print specs
    #p.loglog(ks, 1000*specs,lw=2,marker='o',color=colors[ind])

#p.xlim(.04,10)
p.xlim(.019,7)
p.ylim(5e6,1e11)
p.ylabel(r'$\Delta^2(k)\ [{\rm mK}^2]$')
p.xlabel(r'$k_{\parallel}\ [h{\rm Mpc}^{-1}]$')
params = {'legend.fontsize' : 12}
p.rcParams.update(params)
p.legend()
p.show()

