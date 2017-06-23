#! /usr/bin/env python
import numpy as np,sys
from matplotlib.pyplot import *

threshold = 0.975

F = np.load(sys.argv[1])
#load the bootstrap file
pCvs_boot = F['pCvs']
pCs_boot = F['pCs']
kpl = F['kpl']
freqs = F['freqs']
cmd = F['cmd']

#find the upper limit of the data
pCv_upper = np.percentile(pCvs_boot,97.5,axis=0)

#compute the number of points above the data upper limit
N_above = np.sum(pCs_boot>pCv_upper,axis=0)

err = []
#find the limit for each k
for i in xrange(21):
    j = np.where(N_above[i,:]/np.float(pCvs_boot.shape[0])>threshold)[0].min()
    err.append(pI[i,j])
err = np.array(pk)

#repeat to compute the folded limits
#cut it down the middle and stack it bootwise?
pCvs_boot = np.vstack([pCvs_boot[:10],pCvs_boot[10:]])
pCs_boot = np.vstack([pCs_boot[:10],pCvs_boot[10:]]

#find the upper limit of the data
pCv_upper = np.percentile(pCvs_boot,97.5,axis=0)

#compute the number of points above the data upper limit
N_above = np.sum(pCs_boot>pCv_upper,axis=0)

err_fold = []
#find the limit for each k
for i in xrange(21):
    j = np.where(N_above[i,:]/np.float(pCvs_boot.shape[0])>threshold)[0].min()
    err_fold.append(pI[i,j])
err_fold = np.array(err_fold)
pk_fold = np.ones_like(err_fold)*1e-3


n.savez('pspec.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd,freq=freq)
