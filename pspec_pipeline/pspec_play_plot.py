#! /usr/bin/env python
import numpy as n, pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
import glob
import optparse, sys

#o = optparse.OptionParser()
#o.add_option('--output', type='string', default=None,
#    help='Output filename to save power spectrum values.')
#opts,args = o.parse_args(sys.argv[1:])

count=1
labels = []
labels.append("")
dirs = []
for dir in glob.glob('C_*'):
    if dir in ['C_nofrf','C_sep12','C_sep02','C_original']: 
    #if dir in ['C_original', 'C_reg020', 'C_reg050', 'C_reg100']: #just regularizations
        dirs.append(dir) 
dirs.sort()
for dir in dirs:
    #Read file
    npz = dir+'/factor.npz'
    try: file = n.load(npz)
    except: continue
    print 'Reading', npz.split('/')[0]
    npzI = n.load(dir+'/pspec_I.npz')
    npzC = n.load(dir+'/pspec_C.npz')
    kpl = file['kpl']
    pIvs = npzI['pk']
    pCvs = npzC['pk'] #PS values come from after pspec_cov_boot_v002
    pIvs_err = npzI['err']
    pCvs_err = npzC['err']
    pIv = n.median(n.abs(pIvs))
    pCv = n.median(n.abs(pCvs))
    factor = file['factor']
    #Subplots
    plot1 = p.figure(1,figsize=(8,8))
    numplots = n.ceil(len(dirs)/2.)
    ax = p.subplot(numplots,numplots,count)
    p.errorbar(kpl,pIvs,yerr=2*pIvs_err,fmt='bs',capsize=0,linewidth=1.5,label='pI' if count==1 else "")
    p.errorbar(kpl,n.abs(n.real(pCvs)),yerr=2*pCvs_err,fmt='go',capsize=0,linewidth=1.5,label='pC' if count==1 else "")
    p.errorbar(kpl,n.abs(n.real(pCvs))*factor,yerr=2*pCvs_err*factor,fmt='r^',capsize=0,linewidth=1.5,label='pC corrected' if count==1 else "")
    p.legend(prop={'size':8},loc=4)
    p.xlabel('k')
    p.ylabel('P(k)')
    ax.set_yscale('log',nonposy='clip')
    p.ylim(1e-2,1e14)
    label = npz.split('/')[0].split('_')[1].split('.')[0]
    labels.append(label)
    p.title(label)
    """
    #Overall plot
    p.figure(2)
    p.semilogy(count,pIv,'bs',label='pI' if count==1 else "")
    p.semilogy(count,pCv,'go',label='pC' if count==1 else "")
    p.semilogy(count,pCv*factor,'r^',label='pC corrected' if count==1 else "")
    """
    count += 1
p.figure(1)
p.tight_layout()
"""
p.figure(2)
p.xlim(0,count+1)
p.ylim(1e0,1e14)
xints = n.arange(count+1)
p.xticks(xints,labels)
p.ylabel('Median P(k)')
p.legend(prop={'size':8})
p.title('Playing with C')
"""
p.show()
