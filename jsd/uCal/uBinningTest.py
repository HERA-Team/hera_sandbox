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

#Note to self: I'm not sure If I've implemented the Bayesian information criterion correctly, but in general I think there's no strong evidence for picking anything other than delta u = .5

#############################################
#   Set-up script specific to PAPER 128. 
#   TODO: Hardcoded for now.
#############################################


regenerateEverything = False

dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.57058.xx.uvcRRE.npz")] 
pol='xx'
alwaysFlaggedChannels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 75, 76, 77, 169, 183, 185, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202]
alsoFlagTheseChannels = [14, 15, 55, 101, 123, 179, 180, 181, 182, 184, 186]    
flaggedChannels = sorted(list(set(alwaysFlaggedChannels).union(set(alsoFlagTheseChannels))))
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

uBinSizes = np.arange(.25,3,.25)
results = [{} for i in range(len(uBinSizes))]

for n,uSize in enumerate(uBinSizes):
    print 'Now working on uBinSize = ', uSize
    if regenerateEverything: #regenerate uReds and data
        uReds = uc.uCalReds(freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau=.3) #just pass in freqs
        uReds.applyuCut(uMin=25, uMax=150)
        uReds.applyChannelFlagCut(flaggedChannels) 
        uCal = uc.uCalibrator(uReds.getBlChanPairs())
        uCal.computeVisibilityCorrelations(data, samples)
        pickle.dump([uReds, uCal], open('./Data/uCalData.p', 'wb'))
    
    uReds, uCal = pickle.load(open('./Data/uCalData.p','rb'))
    uCal.setupBinning(uBinSize = uSize, duBinSize = 5.0/203)
    
    print 'Now performing logcal...'
    betasLogcal, SigmasLogcal, DsLogcal = uCal.performLogcal()
    noiseCovDiag = uCal.generateNoiseCovariance(betasLogcal)
    betas, Sigmas, Ds = betasLogcal.copy(), SigmasLogcal.copy(), DsLogcal.copy()
    
    print 'Now performing lincal...'
    previousChiSqPerDoF = 1e10
    for iteration in range(50): 
        betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
        print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
        if np.abs(previousChiSqPerDoF - chiSqPerDoF)/chiSqPerDoF < 1e-4: break
        previousChiSqPerDoF = chiSqPerDoF
        #noiseCovDiag = uCal.generateNoiseCovariance(betas) #updated each cycle based on improved result for beta
    
    #noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
    #betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)    
    print 'Final, noise-median-renormalized chi^2/dof = ' + str(chiSqPerDoF)

    results[n]['chiSqPerDoF'] = chiSqPerDoF
    results[n]['betas'] = betas
    results[n]['bandpassErrors'] = freqs[uCal.chans]**2.55 * np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    results[n]['bandpass'] = np.abs(freqs[uCal.chans]**2.55 * betas)
    results[n]['noiseCovDiag'] = noiseCovDiag

    errors = uCal.computeErrors(betas, Sigmas, Ds)
    nParam = 2*(len(betas)+len(Sigmas)+len(Ds))
    logLike = -.5 * np.sum(np.abs(errors)**2/noiseCovDiag) + -.5 * uCal.nPairs *np.log(np.pi*2)  + -.5 * np.sum(np.log(noiseCovDiag))
    results[n]['BIC'] = -2*logLike + nParam * np.log(uCal.nPairs)
    
#%% 

plt.figure(1); plt.clf()

standardPoly = np.polyfit(np.arange(.1,.2,.1/203)[uCal.chans], np.mean(np.asarray([result['bandpass'] for result in results]),axis=0), 2)
standard = np.polyval(standardPoly, np.arange(.1,.2,.1/203)[uCal.chans])

for result in results[0:6]:
    plt.errorbar(np.arange(.1,.2,.1/203)[uCal.chans],result['bandpass']/standard, result['bandpassErrors']/standard)

plt.xlabel('Frequency (GHz)'); plt.ylabel('Bandpass')
plt.legend(uBinSizes[0:6])

#%%
plt.figure(2); plt.clf()
plt.plot(uBinSizes, [result['chiSqPerDoF'] for result in results])
plt.xlabel('uBinSize'); plt.ylabel('chiSqPerDoF')

#%%
plt.figure(3); plt.clf()
for result in results:
    plt.plot(result['bandpass'])
plt.legend(uBinSizes)
plt.legend(uBinSizes)