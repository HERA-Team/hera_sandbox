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

#############################################
#   Set-up script specific to PAPER 128. 
#   TODO: Hardcoded for now.
#############################################


regenerateEverything = False

dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.43835.xx.uvcRRE.npz")] 
pol='xx'
#alwaysFlaggedChannels = []#[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 75, 76, 77, 169, 183, 185, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202]
flaggedChannels = []#sorted([14, 15, 55, 101, 123, 179, 180, 181, 182, 184, 186,     20,78,170])
freqs = np.arange(.1,.2,.1/203)

if regenerateEverything:
    seps = np.arange(1,16) #* 15m  = baseline lengths
    dx = 15.0 / scipy.constants.c * 1e9
    separations = dx*seps
    #freqs = np.arange(.1,.2,.1/203)[160:]
    chans = range(len(freqs))
    chan2FreqDict = {chan: freqs[chan] for chan in chans}
    uTolerance = 15.0 * 15 * (freqs[1]-freqs[0]) / 3e8 * 1e9  

    aa = a.cal.get_aa('psa6622_v003', freqs)
    blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
    m, _, mdl,_= omni.from_npz(dataFiles[0])
    bls = [(0,i) for i in seps]
    for ij in mdl['xx'].keys():
        bl = a.miriad.ij2bl(ij[0],ij[1]) 
        try:
            s = tuple(map(int,bl2sep[bl].split(',')))
        except(KeyError):
            continue
        if s in bls and s[0] == 0:
            bls[s[-1]-1] = ij
    bl2SepDict = {bl: np.asarray([sep,0.0]) for bl,sep in zip(bls,separations)}
    #bl2SepDict = {bl: np.asarray([sep]) for bl,sep in zip(bls,separations)}

    def blRedundancyDictFromCalFile(calFile = 'psa6622_v003'):
        """ Uses Omnical to figure out the number of redundant baselines for a given separation (as demarkated by a representative baseline tuple). 
        Returns redundancyDict and also saves it as an instance variable."""
        aa = a.cal.get_aa(calFile,np.array([.15]))
        print 'Now generating baseline redundancy dictionary from calfile...'
        ex_ants = [5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110]
        info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants)
        reds = info.get_reds()
        redundancyDict = {}
        #This is overkill...it maps all baselines to their redundancy, not just the ones we need. But that's fine.
        for red in reds:
            for bl in red: redundancyDict[bl] = len(red)
        return redundancyDict

    redundancyDict = blRedundancyDictFromCalFile(calFile = 'psa6622_v003')

    def loadVisibilitiesAndSamples(dataFiles, pol, blList, redundancyDict):
        """Based on Zaki and Carina's omni_get_k.py"""
        print 'Now reading all data files...'
        data = {}
        samples = {}
        jds = []
        lsts = []
        files = []
        conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]

        for fl in dataFiles:
            print '   Reading %s'%fl
            meta,_,mdl,_ = omni.from_npz(fl)
            jd = meta['jds']
            lst = meta['lsts']
            jds.append(jd)
            lsts.append(lst)
            files.append(fl.split('/')[-1]) #ignore the folder
            for b in bls:
                if b in conjugate: _d = np.conj(mdl[pol][b])
                else: _d = mdl[pol][b]
                if b not in data.keys(): data[b] = _d
                else: data[b] = np.concatenate((data[b],_d), axis=0)
                samples[b] = redundancyDict[b] * np.ones(data[b].shape)
                samples[b][data[b]==0] = 0 #TODO: we need a more sophisticated way of handling flagging than this!!!
        return data, samples

    data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict)
        


#############################################
#   uCal Script
#############################################

print '\nNow performing uCal...\n'

if regenerateEverything: #regenerate uReds and data
    uReds = uc.uCalReds(freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau=.3, verbose=True) #just pass in freqs
    uReds.applyuCut(uMin=25, uMax=150)
    uReds.applyChannelFlagCut(flaggedChannels) 
    uCal = uc.uCalibrator(uReds.getBlChanPairs())
    uCal.computeVisibilityCorrelations(data, samples, verbose = True)
    pickle.dump([uReds, uCal], open('./Data/uCalData.p', 'wb'))
else:
    uReds, uCal = pickle.load(open('./Data/uCalData.p','rb'))

while True:
    uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)

    print 'Now performing logcal...'
    betasLogcal, SigmasLogcal, DsLogcal = uCal.performLogcal()
    noiseCovDiag = uCal.generateNoiseCovariance(betasLogcal)
#    noiseCovDiag = np.ones(noiseCovDiag.shape)s
    betas, Sigmas, Ds = betasLogcal.copy(), SigmasLogcal.copy(), DsLogcal.copy()

    print 'Now performing lincal...'
    prevBetas, prevSigmas, prevDs = 0, 0, 0
    for iteration in range(50): 
        betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
        print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
        if np.average(np.abs(betas - prevBetas)/np.abs(betas)) < 1e-4: break
        prevBetas, prevSigmas, prevDs = betas, Sigmas, Ds 
        noiseCovDiag = uCal.generateNoiseCovariance(betas) #updated each cycle based on improved result for beta
#        noiseCovDiag = np.ones(noiseCovDiag.shape)

    noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
    badChans = uCal.identifyBadChannels(betas, Sigmas, Ds, noiseCovDiag, maxAvgError = 2.5)
    if len(badChans)==0: break
    print 'Removing bad channels: ', badChans
    uCal.applyChannelFlagCut(badChans)


betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)    
print 'Final, noise-median-renormalized chi^2/dof = ' + str(chiSqPerDoF)

#%%
def SortedEigensystem(matrix):
    """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
    evals,evecs = np.linalg.eig(matrix)
    indices = np.argsort(np.abs(evals))[::-1]   
    return evals[indices], evecs[:,indices]    

#%% Looking at lincal
params = uCal.nChans + uCal.nuBins + uCal.nduBins
evals, evecs = SortedEigensystem(uCal.AtNinvA)
nullspace = np.asarray([np.append(evecs[0::2,n],evecs[1::2,n]) for n in range(-4,0)])
#fullNullspace = (nullspace[0:1,:]+nullspace[1:2,:]+nullspace[2:3,:]+nullspace[3:4,:]).dot(nullspace[0:1,:]+nullspace[1:2,:]+nullspace[2:3,:]+nullspace[3:4,:]).T.conj())
fullNullspace = np.real(nullspace.conj().T.dot(nullspace))


plt.figure(1); plt.clf()
plt.imshow(fullNullspace)
plt.colorbar()

#%% looking at logcal
params = uCal.nChans + uCal.nuBins + uCal.nduBins
evals, evecs = SortedEigensystem(uCal.logcalBtB)
nullspace = np.asarray([evecs[:,n] for n in range(-2,0)])

fullNullspace = np.real(nullspace.conj().T.dot(nullspace))
v = np.asarray([1]*uCal.nChans + [0]*(uCal.nuBins + uCal.nduBins)) / (uCal.nChans)**.5
v.shape = (params,1)
Proj = np.eye(params) - v.dot(v.T)
newNullspace = (Proj.dot(fullNullspace)).dot(Proj)
evals2, evecs2 = SortedEigensystem(newNullspace)
finalEvec = evecs2[:,0]

plt.figure(2); plt.clf()
plt.plot(finalEvec)
#plt.colorbar()

#note to self: the conclusion of this is that the fourth degeneracy is multiplying all Sigmas by e^i*phi and all Ds by e^i*-phi where phi is constant


