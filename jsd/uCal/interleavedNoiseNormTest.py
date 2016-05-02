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


def runAll(dataFiles, flaggedChannels, findFlaggedChannels=False, subset = None):

    pol='xx'
    freqs = np.arange(.1,.2,.1/203)

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
                if subset is not None:
                    if b not in data.keys(): data[b] = _d[subset::2,:]
                    else: data[b] = np.concatenate((data[b],_d[subset::2,:]), axis=0)
                else:
                    if b not in data.keys(): data[b] = _d
                    else: data[b] = np.concatenate((data[b],_d), axis=0)                    
                samples[b] = redundancyDict[b] * np.ones(data[b].shape)
                samples[b][data[b]==0] = 0 #TODO: we need a more sophisticated way of handling flagging than this!!!
        return data, samples

    data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict)

    print '\nNow performing uCal...\n'


    uReds = uc.uCalReds(freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau=.3, verbose=True) #just pass in freqs
#    uReds.applyuCut(uMin=25, uMax=150)
    uReds.applyChannelFlagCut(flaggedChannels) 
    uCal = uc.uCalibrator(uReds.getBlChanPairs())
    uCal.computeVisibilityCorrelations(data, samples, verbose = True)

    while True:
        uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)

        print 'Now performing logcal...'
        betasLogcal, SigmasLogcal, DsLogcal = uCal.performLogcal()
        noiseCovDiag = uCal.generateNoiseCovariance(betasLogcal)
        betas, Sigmas, Ds = betasLogcal.copy(), SigmasLogcal.copy(), DsLogcal.copy()

        print 'Now performing lincal...'
        prevBetas, prevSigmas, prevDs = 0, 0, 0
        for iteration in range(50): 
            betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
            print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
            if np.average(np.abs(betas - prevBetas)/np.abs(betas)) < 1e-4: break
            prevBetas, prevSigmas, prevDs = betas, Sigmas, Ds 
            noiseCovDiag = uCal.generateNoiseCovariance(betas) #updated each cycle based on improved result for beta

        noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
        badChans = uCal.identifyBadChannels(betas, Sigmas, Ds, noiseCovDiag, maxAvgError = 2.5)        
        if len(badChans)==0 or not findFlaggedChannels: break
        print 'Removing bad channels: ', badChans
        uCal.applyChannelFlagCut(badChans)





    betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)    
    print 'Final, noise-median-renormalized chi^2/dof = ' + str(chiSqPerDoF)

        #Save Results
    print 'Saving results to', dataFiles[-1].replace('.npz','.uCal.npz')
    bandpass, weights = np.zeros(len(freqs),dtype=complex), np.zeros(len(freqs))
    bandpass[uCal.chans] = freqs[uCal.chans]**2.55 * betas
    weights[uCal.chans] = 1
    uc.save2npz(dataFiles[-1].replace('.npz','.uCal.npz'), dataFiles, bandpass, weights, betas, uCal.chans, Sigmas, uCal.uBins, Ds, uCal.duBins, noiseCovDiag, uCal.A)

    flaggedChannels = list(set(chans).difference(set(list(uCal.chans))))
    print flaggedChannels
    print str(len(set(flaggedChannels)))
    return betas, Sigmas, Ds, flaggedChannels, uCal, noiseCovDiag




#############################################
#   Now do everything
#############################################



dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
betas, Sigmas, Ds, flaggedChannels, uCal, noiseCovDiag = runAll(dataFiles, [], findFlaggedChannels=True)

betasEven, SigmasEven, DsEven, dummy, uCalEven, dummy = runAll(dataFiles, flaggedChannels, subset = 0)
betasOdd, SigmasOdd, DsOdd, dummy, uCalOdd, dummy = runAll(dataFiles, flaggedChannels, subset = 1)

#%%##########################################
#   Now Plot
#############################################

allChans = sorted(list(set(range(203)).difference(set(flaggedChannels))))

plt.figure(1); plt.clf()
plt.plot(allChans, np.abs(betas),'.')
plt.plot(allChans, np.abs(betasEven),'.')
plt.plot(allChans, np.abs(betasOdd),'.')
#plt.plot(np.abs(Sigmas),'.')
#plt.plot(np.abs(SigmasEven),'.')
#plt.plot(np.abs(SigmasOdd),'.')
plt.legend(['Both','Even','Odd'])
plt.xlabel('abs(beta)')
plt.ylabel('chan')

evenSamples = np.sum(np.asarray([item['samples'] for item in uCalEven.blChanPairs.values()]))
oddSamples = np.sum(np.asarray([item['samples'] for item in uCalOdd.blChanPairs.values()]))
samples = np.sum(np.asarray([item['samples'] for item in uCal.blChanPairs.values()]))
print evenSamples, oddSamples, samples