import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import scipy
import scipy.constants
import omnical
import time
import aipy as a
import optparse, sys, os
import capo
import capo.uCal as uc


allResultsFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.uCal.npz")]
results = [np.load(resultFile) for resultFile in allResultsFiles]

#%%
plt.figure(1); plt.clf()
allBandpasses = np.asarray([result['bandpass'] for result in results])
plt.plot(np.abs(np.average(allBandpasses,axis=0)))

#%%

allChans = set(results[0]['chans'])
alluBins = set([tuple(uBin) for uBin in results[0]['uBins']])
allduBins = set(results[0]['duBins'])
for result in results:
    allChans = allChans.union(result['chans'])
    alluBins = alluBins.union(set([tuple(uBin) for uBin in result['uBins']]))
    allduBins = allduBins.union(result['duBins'])
allChans, alluBins, allduBins = [sorted(list(thisSet)) for thisSet in (allChans, alluBins, allduBins)]
nChans, nuBins, nduBins = len(allChans), len(alluBins), len(allduBins)
nParams = nChans + nuBins + nduBins

#%%
invCovs, xHats = [], []
for result in results:
    chanIndices = [(allChans.index(chan))*2 for chan in result['chans']]
    uBinIndices = [(alluBins.index(tuple(uBin)) + nChans)*2 for uBin in result['uBins']]
    duBinIndices = [(allduBins.index(duBin) + nChans + nuBins)*2 for duBin in result['duBins']]
    allIndices = chanIndices + uBinIndices + duBinIndices
    noiseCovDiag = result['noiseCovDiag']
    A = result['A'][()]
    Ninv = csr_matrix((np.append((noiseCovDiag)**-1,(noiseCovDiag)**-1), (np.arange(2*len(noiseCovDiag)), np.arange(2*len(noiseCovDiag)))))
    invCovHere = ((A.conjugate().T.dot(Ninv)).dot(A)).toarray()
    xHatHere = np.asarray([[np.real(x), np.imag(x)] for x in np.append(result['betas'], np.append(result['Sigmas'], result['Ds']))]).flatten()     
    invCov = np.zeros((nParams*2,nParams*2))
    xHat = np.zeros((nParams*2,))
    for i,ind1 in enumerate(allIndices):
        xHat[ind1] = xHatHere[2*i]
        xHat[ind1+1] = xHatHere[2*i+1]        
        for j,ind2 in enumerate(allIndices):
            invCov[ind1,ind2] = invCovHere[2*i,2*j]
            invCov[ind1+1,ind2] = invCovHere[2*i+1,2*j]
            invCov[ind1,ind2+1] = invCovHere[2*i,2*j+1]
            invCov[ind1+1,ind2+1] = invCovHere[2*i+1,2*j+1]  
    invCovs.append(invCov)
    xHats.append(xHat)

        
#%%
CinvSum, CinvxSum = None, None
for Cinv, xHat in zip(invCovs, xHats):

#    Cinv = np.identity(len(Cinv))
#    for i in range(len(xHat)): 
#        if xHat[i]==0: Cinv[i,i]=0
    
    if CinvSum is None:
        CinvSum = Cinv
        CinvxSum  = Cinv.dot(xHat)
    else:
        CinvSum += Cinv
        CinvxSum += Cinv.dot(xHat)

xHat = np.linalg.pinv(CinvSum).dot(CinvxSum)



#%
combinedBandpass = xHat[0:2*nChans:2] + 1.0j*xHat[1:2*nChans:2]
plt.figure(2); plt.clf()
plt.plot(allChans, np.abs(combinedBandpass),'.')
#plt.plot(results[9]['chans'],np.abs(results[9]['betas']),'.')
#%%
plt.figure(222); plt.clf()
for Cinv, xHat in zip(invCovs, xHats):
    print np.max(np.abs(Cinv.dot(xHat)))
    plt.semilogy(np.abs(Cinv.dot(xHat)))

#%%
plt.figure(223); plt.clf()
for result in results:
   plt.plot(result['chans'],np.abs(result['betas']))