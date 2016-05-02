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


#TODO: look into the effects of the range of u values included
#TODO: add in polynomial fitting of final beta solution
#TODO: look at phase ramp
#TODO: root median sq test

#TODO: take N bootstraps (length all the data but with replacement) and use that as a new noise model in noiseCovDiag
#TODO: stop removing high variance channels

#%%##########################################
#   Set-up functions specific to PAPER 128
#############################################

def getBaselines(freqs, intSeps, exampleDataFile, calFile='psa6622_v003'):
    aa = a.cal.get_aa(calFile, freqs)
    blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
    m, _, mdl,_= omni.from_npz(exampleDataFile)
    bls = [(0,i) for i in intSeps]
    for ij in mdl['xx'].keys():
        bl = a.miriad.ij2bl(ij[0],ij[1]) 
        try:
            s = tuple(map(int,bl2sep[bl].split(',')))
        except(KeyError):
            continue
        if s in bls and s[0] == 0:
            bls[s[-1]-1] = ij
    return bls

def blRedundancyDictFromCalFile(calFile = 'psa6622_v003', verbose = False):
    """ Uses Omnical to figure out the number of redundant baselines for a given separation (as demarkated by a representative baseline tuple). 
    Returns redundancyDict and also saves it as an instance variable."""
    aa = a.cal.get_aa(calFile,np.array([.15]))
    if verbose: print 'Now generating baseline redundancy dictionary from calfile...'
    ex_ants = [5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110]
    info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants)
    reds = info.get_reds()
    redundancyDict = {}
    #This is overkill...it maps all baselines to their redundancy, not just the ones we need. But that's fine.
    for red in reds:
        for bl in red: redundancyDict[bl] = len(red)
    return redundancyDict

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

def splitDataAndSamples(data, samples, bls, nWaySplit):
    """Splits data and samples n ways, alternating. Returns a list of data dictionaries and a list of samples dictionaries."""
    splitData, splitSamples = [{} for i in range(nWaySplit)], [{} for i in range(nWaySplit)]
    for bl in bls:
        for i in range(nWaySplit):
            splitData[i][bl] = data[bl][i::nWaySplit,:]
            splitSamples[i][bl] = samples[bl][i::nWaySplit,:]
    return splitData, splitSamples


def harmonizeChannelFlags(chans, uCal, uCalSplits, verbose=False):
    """This function flags all channels that are flagged in any of the splits for both uCal and all splits."""
    unflagged = set(uCal.chans)
    for split in uCalSplits: unflagged &= set(split.chans)
    toFlag = list(set(chans) - unflagged)
    uCal.applyChannelFlagCut(toFlag)
    for split in uCalSplits: split.applyChannelFlagCut(toFlag)
    if verbose: print "After harmonizing channels flags, " + str(len(toFlag)) + " channels have been flagged for all split parts."

def computeSplitVarianceNormalization(betas, splitBetas, uCal, uCalSplits):
    """This function uses the observed variance between the splits to normalize the noise covariance."""
    totalSamples = np.sum([value['samples'] for value in uCal.blChanPairs.values()])
    convertVarFactor = totalSamples * np.sum([1.0/np.sum([value['samples'] for value in uCalSplits[i].blChanPairs.values()]) for i in range(nWaySplit)]) / len(uCalSplits)
    overallBandpassVariance = np.mean(np.abs(betas)**2 * ((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2])))
    overallObservedVariance = np.mean(np.var(np.asarray(splitBetas),axis=0)/convertVarFactor)
    return overallObservedVariance/overallBandpassVariance
    
#############################################
#   uCal Parameters
#############################################

regenerateEverything = True
verbose = True
dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.43835.xx.uvcRRE.npz")] 
nWaySplit = 19

#%%##########################################
#   uCal Setup
#############################################

#PAPER 128 Setup
pol='xx'
freqs = np.arange(.1,.2,.1/203)
chans = range(len(freqs))
chan2FreqDict = {chan: freqs[chan] for chan in chans}
intSeps = np.arange(1,16) #* 15m  = baseline lengths
dx = 15.0 / scipy.constants.c * 1e9
separations = dx*intSeps
bls = getBaselines(freqs, intSeps, dataFiles[0], calFile='psa6622_v003')
bl2SepDict = {bl: np.asarray([sep,0.0]) for bl,sep in zip(bls,separations)}
redundancyDict = blRedundancyDictFromCalFile(calFile = 'psa6622_v003')

if regenerateEverything:
    uReds = uc.uCalReds(freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau=.3, verbose=verbose) #just pass in freqs
    uReds.applyuCut(uMin=25, uMax=150)
    #uReds.applyChannelFlagCut(flaggedChannels) 
    uCal = uc.uCalibrator(uReds.getBlChanPairs())
    data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict)
    uCal.computeVisibilityCorrelations(data, samples, verbose=verbose)
    uCalSplits = [uc.uCalibrator(uReds.getBlChanPairs()) for i in range(nWaySplit)]
    splitData, splitSamples = splitDataAndSamples(data, samples, bls, nWaySplit)
    for i in range(nWaySplit): uCalSplits[i].computeVisibilityCorrelations(splitData[i], splitSamples[i], verbose=verbose)
    harmonizeChannelFlags(chans, uCal, uCalSplits, verbose=verbose)
    pickle.dump([uReds, uCal, uCalSplits, dataFiles], open('./Data/uCalDataWithSplit.p', 'wb'))
else: 
    uReds, uCal, uCalSplits, loadedDataFiles = pickle.load(open('./Data/uCalDataWithSplit.p','rb'))
    if not len(uCalSplits) == nWaySplit or not loadedDataFiles == dataFiles: raise RuntimeError('Loaded uCalData.p does not match current specifications.')

#%%##########################################
#   Run uCal
#############################################

def performuCal(uCal, maxIterations = 20, verbose = False, initialGuesses = None):
    """Performs logcal and then up to maxIterations of lincal, returning betas, Sigmas, Ds, and the noise covariance diagonal. If inital gueses are provided, logcal is skipped."""
    uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)
    if verbose: print 'Now performing logcal...'
    betasLogcal, SigmasLogcal, DsLogcal = uCal.performLogcal()
    betas, Sigmas, Ds = betasLogcal.copy(), SigmasLogcal.copy(), DsLogcal.copy()
    if initialGuesses is not None: betas, Sigmas, Ds = initialGuesses
    noiseCovDiag = uCal.generateNoiseCovariance(betas)

    if verbose: print 'Now performing lincal...'
    prevBetas, prevSigmas, prevDs = 0, 0, 0
    for iteration in range(maxIterations): 
        betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
        if verbose: print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
        if np.average(np.abs(betas - prevBetas)/np.abs(betas)) < 1e-4: break
        prevBetas, prevSigmas, prevDs = betas, Sigmas, Ds 
        noiseCovDiag = uCal.generateNoiseCovariance(betas) #updated each cycle based on improved result for beta
    noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
    return betas, Sigmas, Ds, noiseCovDiag

#%%

while True:
    betas, Sigmas, Ds, noiseCovDiag = performuCal(uCal, verbose=verbose)
    badChans1 = uCal.identifyBadChannels(betas, Sigmas, Ds, noiseCovDiag, maxAvgError = 10)
    uCal.applyChannelFlagCut(badChans1)
    harmonizeChannelFlags(chans, uCal, uCalSplits, verbose=verbose)
    if verbose: print 'Removed ' + str(len(badChans1)) + ' bad channels for bad fitting:', badChans1

    betas, Sigmas, Ds, noiseCovDiag = performuCal(uCal, verbose=verbose, initialGuesses=[betas, Sigmas, Ds])
    splitBetas, splitSigmas, splitDs, splitNoiseCovDiags = [[None for i in range(nWaySplit)] for j in range(4)]    
    for i in range(nWaySplit):
        splitBetas[i], splitSigmas[i], splitDs[i], splitNoiseCovDiags[i] = performuCal(uCalSplits[i], verbose=verbose, initialGuesses=[betas, Sigmas, Ds])
    modelVarianceByChannel = np.abs(betas)**2 * ((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))
    observedVarianceByChannel = np.var(np.asarray(splitBetas),axis=0)/nWaySplit
    noiseCovDiag *= np.median(observedVarianceByChannel) / np.median(modelVarianceByChannel)
    betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)        
    
    badChans2 = np.asarray(uCal.chans)[(observedVarianceByChannel/modelVarianceByChannel) > 10]
    uCal.applyChannelFlagCut(badChans2)
    harmonizeChannelFlags(chans, uCal, uCalSplits, verbose=verbose)
    if verbose: print 'Removed ' + str(len(badChans2)) + ' bad channels for too high observed variance:', badChans2
    if len(badChans1)==0 and len(badChans2)==0: break

print 'Final, noise-median-renormalized chi^2/dof = ' + str(chiSqPerDoF)


#%%


plt.figure(11); plt.clf()
modelVariance = np.abs(betas)**2 * ((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))
plt.errorbar(uCal.chans, modelVariance,yerr=modelVariance * ((nWaySplit-3.0)/nWaySplit/(nWaySplit-1))**.5)
plt.semilogy(uCal.chans, np.var(np.asarray(splitBetas),axis=0)/nWaySplit,'.')
plt.legend(['Observed Variance','Modeled variance'])
plt.xlabel('chan')

plt.figure(12); plt.clf()
modeledVariance = np.abs(betas)**2 * ((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))
observedVariance = np.var(np.asarray(splitBetas),axis=0)/nWaySplit
plt.plot(uCal.chans, observedVariance/modeledVariance)
plt.xlabel('chan')
plt.ylabel('observed variance')

#%%
#
#"""TODO: 
#Take the median value of the ratio of the observed variance to the noise.
#Use that to renormalize the noise model
#Now look for outliers relative to the noise model
#"""
#
#
##%%
##
#
##%%
#plt.figure(10); plt.clf(); 
#plt.plot(uCal.chans, np.abs(np.asarray(splitBetas).T),'.')
#plt.xlabel('chan'); plt.ylabel('abs(beta)')
#plt.plot(uCal.chans, np.abs(betas))
#
#plt.figure(15); plt.clf()
#splitRelativeStd = np.std(np.asarray(splitBetas),axis=0) / np.abs(betas)
#plt.plot(uCal.chans, splitRelativeStd,'.')
#plt.xlabel('chan')
#plt.title('std(interleaved betas) / abs(full beta)')



#############################################
#   Diagnostic Plotting
#############################################
if True:
     #%% Setup
     def lincalScatter(x,y,color=None, figNum=100, xs='log', ys='log', title='', clear=True, xl='', yl=''):
         plt.figure(figNum); 
         if clear: plt.clf()
         if color is not None: plt.scatter(x,y,c=color)
         else: plt.scatter(x,y)
         plt.yscale(ys); plt.xscale(xs)
         plt.ylim([.9*np.min(y), 1.1*np.max(y)]); plt.xlim([.9*np.min(x), 1.1*np.max(x)])
         plt.xlabel(xl); plt.ylabel(yl); plt.title(title)
     duList = np.asarray([entry['du'] for entry in uCal.blChanPairs.values()])
     uList = np.asarray([entry['u'] for entry in uCal.blChanPairs.values()])
     ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
     ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
     relativeScatter = np.zeros(len(freqs))
     relativeScatter[uCal.chans] = observedVarianceByChannel/modelVarianceByChannel
     relativeScatterList = np.asarray([relativeScatter[ch1] + relativeScatter[ch2] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])

    #%%Betas
     plt.figure(101); plt.clf()
     inferredErrorsOnAbsBeta = np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
     plt.errorbar(uCal.chans,np.abs(betas),yerr=inferredErrorsOnAbsBeta)
     plt.plot(uCal.chans, np.abs(np.asarray(splitBetas).T),':.')#'--.')
     plt.xlabel('chan')#     plt.xlabel('Frequency (GHz)'); 
     plt.ylabel('Abs(Beta)'); plt.title('Beta vs. Split Observation Beta')

     #%%Bandpass
     plt.figure(1); plt.clf()
     inferredErrorsOnAbsBeta = freqs[uCal.chans]**2.55 * np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
     plt.errorbar(np.arange(.1,.2,.1/203)[uCal.chans],np.abs(freqs[uCal.chans]**2.55 * betas),yerr=inferredErrorsOnAbsBeta)
#     plt.errorbar(uCal.chans,np.abs(freqs[uCal.chans]**2.55 * betas),yerr=inferredErrorsOnAbsBeta)
#     plt.plot(np.arange(.1,.2,.1/203)[uCal.chans], np.asarray([freqs[uCal.chans]**2.55 * beta for beta in np.abs(np.asarray(splitBetas))]).T,'.')
     plt.xlabel('Frequency (GHz)'); 
     plt.ylabel('Abs(Lincal Bandpass)'); plt.title(r'Lincal Bandpass Corrected by $\nu^{2.55}$')

     #%%Predicted vs. Observed Scatter
     plt.figure(2); plt.clf()
     visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
     predictedCorrs = visCorrs - uCal.computeErrors(betas, Sigmas, Ds)
     lincalScatter(np.abs(predictedCorrs), np.abs(visCorrs), color=np.linalg.norm(uList,axis=1), figNum=2, ys='log', xs='log', title = 'color = u') #color=np.log10(relativeScatterList)
     plt.plot([0,1],[0,1],'k--')
     plt.xlabel('Abs(Predicted Correlations)'); plt.ylabel('Abs(Observe Correlations)')

     #%%Examine error correlations
     plt.figure(3); plt.clf()
     AtNinvAinv = np.linalg.pinv(uCal.AtNinvA)[0:2*uCal.nChans:2,0:2*uCal.nChans:2]
     inverseSqrtDiag = np.diag(np.diag(AtNinvAinv)**-.5)
     plt.imshow(inverseSqrtDiag.dot(AtNinvAinv.dot(inverseSqrtDiag)), interpolation='nearest', extent = [uCal.chans[0],uCal.chans[-1],uCal.chans[0],uCal.chans[-1]], vmin=-.2, vmax=1)
     plt.title('Error Correlation Matrix for Real Part of Beta')
     plt.xlabel('Channel'); plt.ylabel('Channel')
     plt.colorbar()

     #%%Identify bad channels
     chanCompiledList = {chan: [] for chan in uCal.chans}
     errorList = uCal.computeErrors(betas, Sigmas, Ds)
     for f1,f2,error,Nii in zip(ch1List,ch2List,errorList,noiseCovDiag):
            chanCompiledList[f1].append(error/((2*Nii)**.5))
            chanCompiledList[f2].append(error/((2*Nii)**.5))
     chanAvgErrors = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList[chan]))**2)**.5 for chan in uCal.chans])
     plt.figure(4); plt.clf()
     plt.plot(uCal.chans, chanAvgErrors,'.')
     plt.ylabel('Channel-Averaged, Noise Weighted Errors'); plt.xlabel('Channel')
#     badChans = np.asarray(uCal.chans)[chanAvgErrors > 2.5]
#     if len(badChans) > 0: print 'Channels with average sigma > 2.5: ', badChans

print 'Flagged channels: ', len(freqs) - uCal.nChans


