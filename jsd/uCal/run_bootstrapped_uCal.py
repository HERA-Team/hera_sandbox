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
from joblib import Parallel, delayed
import multiprocessing
from scipy.optimize import curve_fit

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

def dataAndSamplesBootstraps(data, samples, bls, nBootstraps):
    """This function returns a list of """
    dataBoostraps, samplesBoostraps = [{} for i in range(nBootstraps)], [{} for i in range(nBootstraps)] 
    nIntegrations = len(data[bls[0]])
    for i in range(nBootstraps):
        indices = np.random.random_integers(0,nIntegrations-1,nIntegrations)        
        for bl in bls:
            dataBoostraps[i][bl] = data[bl][indices,:]
            samplesBoostraps[i][bl] = samples[bl][indices,:]
    return dataBoostraps, samplesBoostraps    

def harmonizeChannelFlags(chans, uCal, uCalSplits, verbose=False):
    """This function flags all channels that are flagged in any of the splits for both uCal and all splits."""
    unflagged = set(uCal.chans)
    for split in uCalSplits: unflagged &= set(split.chans)
    toFlag = list(set(chans) - unflagged)
    uCal.applyChannelFlagCut(toFlag)
    for split in uCalSplits: split.applyChannelFlagCut(toFlag)
    if verbose: print "After harmonizing channels flags, " + str(len(toFlag)) + " channels have been flagged for all bootstraps."

#def computeSplitVarianceNormalization(betas, splitBetas, uCal, uCalSplits):
#    """This function uses the observed variance between the splits to normalize the noise covariance."""
#    totalSamples = np.sum([value['samples'] for value in uCal.blChanPairs.values()])
#    convertVarFactor = totalSamples * np.sum([1.0/np.sum([value['samples'] for value in uCalSplits[i].blChanPairs.values()]) for i in range(nWaySplit)]) / len(uCalSplits)
#    overallBandpassVariance = np.mean(np.abs(betas)**2 * ((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2])))
#    overallObservedVariance = np.mean(np.var(np.asarray(splitBetas),axis=0)/convertVarFactor)
#    return overallObservedVariance/overallBandpassVariance


#############################################
#   uCal Parameters
#############################################

regenerateEverything = False
verbose = True
dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.43835.xx.uvcRRE.npz")] 
nBootstraps = 40
uMin, uMax = 25, 100

#TODO: make this automatic
manualChannelFlags = [14, 17, 55, 74, 93, 101, 102, 152, 153, 166, 168, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186]

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
    #uReds.applyuCut(uMin=25, uMax=150)
    #uReds.applyChannelFlagCut(flaggedChannels) 
    uCal = uc.uCalibrator(uReds.getBlChanPairs())
    data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict)
    uCal.computeVisibilityCorrelations(data, samples, verbose=verbose)
    uCalBoostraps = [uc.uCalibrator(uReds.getBlChanPairs()) for i in range(nBootstraps)]
    dataBoostraps, samplesBoostraps = dataAndSamplesBootstraps(data, samples, bls, nBootstraps)
    for i in range(nBootstraps): uCalBoostraps[i].computeVisibilityCorrelations(dataBoostraps[i], samplesBoostraps[i], verbose=verbose)
    harmonizeChannelFlags(chans, uCal, uCalBoostraps, verbose=verbose)
    pickle.dump([uReds, uCal, uCalBoostraps, dataFiles], open('./Data/uCalDataWithBootstraps.p', 'wb'))
else: 
    uReds, uCal, uCalBoostraps, loadedDataFiles = pickle.load(open('./Data/uCalDataWithBootstraps.p','rb'))
    if not len(uCalBoostraps) == nBootstraps or not loadedDataFiles == dataFiles: raise RuntimeError('Loaded uCalDataWithBootstraps.p does not match current specifications.')

#%%##########################################
#   Run uCal
#############################################

uCal.applyChannelFlagCut(manualChannelFlags)
uCal.applyuCut(uMin=uMin, uMax=uMax)
uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)
for bootstrap in uCalBoostraps: 
    bootstrap.applyChannelFlagCut(manualChannelFlags)    
    bootstrap.applyuCut(uMin=uMin, uMax=uMax)
    bootstrap.setupBinning(uBinSize = .5, duBinSize = 5.0/203)


def performuCal(uCal, maxIterations = 40, verbose = False):
    if not uCal.binningIsSetup: uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)
    if verbose: print 'Now performing logcal...'
    betas, Sigmas, Ds = uCal.performLogcal()
    noiseCovDiag = uCal.generateNoiseCovariance(betas)
    if verbose: print 'Now performing lincal...'
    prevBetas, prevSigmas, prevDs = 0, 0, 0
    for iteration in range(maxIterations): 
        betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
        if verbose: print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
        if np.average(np.abs(betas - prevBetas)/np.abs(betas)) < 1e-6: break
        prevBetas, prevSigmas, prevDs = betas, Sigmas, Ds 
        noiseCovDiag = uCal.generateNoiseCovariance(betas)
    noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
    return betas, Sigmas, Ds, noiseCovDiag


betas, Sigmas, Ds, noiseCovDiag = performuCal(uCal, verbose=verbose)
bootstrapBetas, bootstrapSigmas, bootstrapDs =  [[None for i in range(nBootstraps)] for j in range(3)]    
parallelBoostrapResults = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(performuCal)(uCalBoostraps[i]) for i in range(nBootstraps))
for i in range(nBootstraps): bootstrapBetas[i], bootstrapSigmas[i], bootstrapDs[i], dummy = parallelBoostrapResults[i]

#%%##########################################
#   Iterative Bootstrapping
#############################################
allResults = []
import copy
allResults.append([copy.copy(betas), copy.copy(Sigmas), copy.copy(Ds), copy.copy(noiseCovDiag), copy.copy(bootstrapBetas), copy.copy(bootstrapSigmas), copy.copy(bootstrapDs)])

def performuCal2(uCal, noiseCovDiag, maxIterations = 40, verbose = False):
    if not uCal.binningIsSetup: uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)    
    if verbose: print 'Now performing logcal...'    
    betas, Sigmas, Ds = uCal.performLogcal()
    if verbose: print 'Now performing lincal...'
    prevBetas, prevSigmas, prevDs = 0, 0, 0
    for iteration in range(40): 
        betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)
        if verbose: print '    ' + str(iteration) + ') chi^2/dof = ' + str(chiSqPerDoF)
        if np.average(np.abs(betas - prevBetas)/np.abs(betas)) < 1e-6: break
        prevBetas, prevSigmas, prevDs = betas, Sigmas, Ds 
    return betas, Sigmas, Ds

for i in range(5):
    if verbose: print 'Now working on iteration', i+2
    betasStd = np.std(np.asarray(bootstrapBetas).T,axis=1)
    SigmasStd = np.std(np.asarray(bootstrapSigmas).T,axis=1)
    DsStd = np.std(np.asarray(bootstrapDs).T,axis=1)
    noiseCovDiag = uCal.modelNoiseVariance(betas, Sigmas, Ds, betasStd, SigmasStd, DsStd)
    betas, Sigmas, Ds = performuCal2(uCal, noiseCovDiag, verbose=verbose)
    noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
    
    data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict)
    uCalBoostraps = [uc.uCalibrator(uReds.getBlChanPairs()) for i in range(nBootstraps)]
    dataBoostraps, samplesBoostraps = dataAndSamplesBootstraps(data, samples, bls, nBootstraps)
    for i in range(nBootstraps): uCalBoostraps[i].computeVisibilityCorrelations(dataBoostraps[i], samplesBoostraps[i], verbose=verbose)
    for bootstrap in uCalBoostraps: 
        bootstrap.applyChannelFlagCut(manualChannelFlags)    
        bootstrap.applyuCut(uMin=uMin, uMax=uMax)
        bootstrap.setupBinning(uBinSize = .5, duBinSize = 5.0/203)
    
    bootstrapBetas, bootstrapSigmas, bootstrapDs =  [[None for i in range(nBootstraps)] for j in range(3)]    
    parallelBoostrapResults = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(performuCal2)(uCalBoostraps[i], noiseCovDiag) for i in range(nBootstraps))
    for i in range(nBootstraps): bootstrapBetas[i], bootstrapSigmas[i], bootstrapDs[i] = parallelBoostrapResults[i]
    allResults.append([copy.copy(betas), copy.copy(Sigmas), copy.copy(Ds), copy.copy(noiseCovDiag), copy.copy(bootstrapBetas), copy.copy(bootstrapSigmas), copy.copy(bootstrapDs)])


#%%
#examine allResults
plt.figure(234); plt.clf()
for result in allResults:
    plt.plot(uCal.chans, np.abs(result[0]))

plt.legend(range(len(allResults)))
plt.xlabel('chan'); plt.ylabel('abs(beta)'); plt.title('Iterative Bootstrapped Noise Variance Convergence')

#%%##########################################
#   Complex Fitting
#############################################



length = .1
def fourier(x, *a):
    ret = a[0]
    for deg in range(1, len(a)):
        ret += a[deg] * np.cos((deg) * np.pi / length * x)
    return ret

def cosineFit(x, data, errors, degree):
    popt, pcov = curve_fit(fourier, x, data, sigma=errors, p0=[1.0] * (degree+1))
    return fourier(x, *popt)
    
def BIC(data, realModel, imagModel, realErrors, imagErrors, DoF):
    chiSq = np.sum(((np.real(data) - realModel)/realErrors)**2) + np.sum(((np.imag(data) - imagModel)/imagErrors)**2)
    return np.sum(np.log(2*np.pi*realErrors**2)) + np.sum(np.log(2*np.pi*imagErrors**2)) + chiSq + DoF*np.log(2*len(data))

def findBestModel(x, data, realErrors, imagErrors, degreeRange):
    BICs = []
    for degree in degreeRange:
        realModel = cosineFit(x, np.real(data), realErrors, degree)
        imagModel = cosineFit(x, np.imag(data), imagErrors, degree)
        BICs.append(BIC(data, realModel, imagModel, realErrors, imagErrors, DoF=2*degree+2))
    bestDegree = (degreeRange[BICs==np.min(BICs)][0])
    bestModel = cosineFit(x, np.real(data), realErrors, bestDegree) + 1.0j * cosineFit(x, np.imag(data), imagErrors, bestDegree)
    bestChiSqPerDoF = (np.sum(((np.real(data) - np.real(bestModel))/realErrors)**2) + np.sum(((np.imag(data) - np.imag(bestModel))/imagErrors)**2))/(len(data)*2)
    return bestModel, BICs, bestDegree, bestChiSqPerDoF

degreeRange = np.arange(60,140)
observedRealErrors = np.std(np.real(np.asarray(bootstrapBetas).T),axis=1)
observedImagErrors = np.std(np.imag(np.asarray(bootstrapBetas).T),axis=1)
modelRealErrors = np.real(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]))**.5
modelImagErrors = np.imag(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
model, BICs, bestDegree, bestChiSqPerDoF = findBestModel(freqs[uCal.chans] - .1, betas*freqs[uCal.chans]**2.55, observedRealErrors*freqs[uCal.chans]**2.55, observedImagErrors*freqs[uCal.chans]**2.55, degreeRange)
print 'best degree', bestDegree
print 'chi^2 per dof of fit =', bestChiSqPerDoF

#%% Plot fits
plt.figure(755); plt.clf()
plt.plot(degreeRange,BICs)
plt.xlabel('Degree'); plt.ylabel('BIC')

plt.figure(766); plt.clf()
plt.errorbar(uCal.chans, np.real(betas)*freqs[uCal.chans]**2.55, yerr=observedRealErrors*freqs[uCal.chans]**2.55)
plt.errorbar(uCal.chans, np.imag(betas*freqs[uCal.chans]**2.55),yerr=observedImagErrors*freqs[uCal.chans]**2.55)
plt.plot(uCal.chans, np.real(model),'r')
plt.plot(uCal.chans, np.imag(model),'k')
plt.xlabel('channel'); plt.ylabel('beta')
plt.legend(['Re(beta) fit', 'Im(beta) fit','Re(beta)','Im(beta)'])

#%%##########################################
#   Diagnostic Plotting
#############################################
if True:

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
#     relativeScatter = np.zeros(len(freqs))
#     relativeScatter[uCal.chans] = observedVarianceByChannel/modelVarianceByChannel
#     relativeScatterList = np.asarray([relativeScatter[ch1] + relativeScatter[ch2] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])

    #%%Betas
    plt.figure(102); plt.clf()
    inferredErrorsOnAbsBeta = np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    plt.errorbar(uCal.chans,np.abs(betas),yerr=inferredErrorsOnAbsBeta)
    plt.plot(uCal.chans, np.abs(np.asarray(bootstrapBetas).T),':.')#'--.')
    plt.xlabel('chan')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Beta)'); plt.title('Beta with Bootstraps')

   #%%Sigmas
    plt.figure(103); plt.clf()
    inferredErrorsOnAbsSigma = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]))**.5
    plt.errorbar(np.asarray(uCal.uBins)[:,0], np.abs(Sigmas),yerr=inferredErrorsOnAbsSigma)
    plt.plot(np.asarray(uCal.uBins)[:,0], np.abs(np.asarray(bootstrapSigmas).T),':.')#'--.')
    plt.xlabel('uBin')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Sigma)'); plt.title('Sigma with Bootstraps')     

   #%%Ds
    plt.figure(104); plt.clf()
    inferredErrorsOnAbsD = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans+2*uCal.nuBins :: 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans+2*uCal.nuBins :: 2]))**.5
#    plt.errorbar(uCal.duBins, np.abs(Ds),yerr=inferredErrorsOnAbsD)
#    plt.plot(uCal.duBins, np.abs(np.asarray(bootstrapDs).T),':.')#'--.')
    plt.errorbar(uCal.duBins, np.abs(allResults[-1][2]), yerr=inferredErrorsOnAbsD)
    plt.plot(uCal.duBins, np.abs(np.asarray(allResults[-1][6]).T),':.')#'--.')
    plt.xlabel('uBin')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Ds)'); plt.title('Ds with Bootstraps')   

    #%%Bandpass
    plt.figure(1); plt.clf()
    inferredErrorsOnBandpass = freqs[uCal.chans]**2.55 * np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    plt.errorbar(np.arange(.1,.2,.1/203)[uCal.chans],np.abs(freqs[uCal.chans]**2.55 * betas),yerr=inferredErrorsOnBandpass)
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
    
    #%%
    plt.figure(4); plt.clf()
    plt.plot(uCal.chans, inferredErrorsOnAbsBeta)
    plt.plot(uCal.chans, np.std(np.abs(np.asarray(bootstrapBetas).T),axis=1),'.-')
    plt.legend(['Inferred Errors','Bootstrap Stds'])
    plt.xlabel('Channel')

    plt.figure(5); plt.clf()
    plt.plot(np.asarray(uCal.uBins)[:,0], inferredErrorsOnAbsSigma)
    plt.plot(np.asarray(uCal.uBins)[:,0], np.std(np.abs(np.asarray(bootstrapSigmas).T),axis=1),'.-')
    plt.legend(['Inferred Errors','Bootstrap Stds'])
    plt.xlabel('uBin')

    plt.figure(6); plt.clf()
    plt.plot(uCal.duBins, inferredErrorsOnAbsD)
    plt.plot(uCal.duBins, np.std(np.abs(np.asarray(bootstrapDs).T),axis=1),'.-')
    plt.legend(['Inferred Errors','Bootstrap Stds'])
    plt.xlabel('duBin')
    
    
    #%% Compare errors to various other quantities
    
    #%% Error vs. Noise
    plt.figure(25); plt.clf()
    lincalScatter(noiseCovDiag**.5, np.abs(uCal.computeErrors(betas, Sigmas, Ds)), figNum=25, color=np.linalg.norm(uList,axis=1))
    plt.plot([0,1],[0,1],'k--')
    plt.xlabel('Expected Std'); plt.ylabel('Abs(Error on VisCorr)'); plt.title('Noise Model Based on Bootstrapped Variances (color = u)')
    
    plt.figure(26); plt.clf()
    lincalScatter(uCal.generateNoiseCovariance(betas)**.5, np.abs(uCal.computeErrors(betas, Sigmas, Ds)), figNum=26, color=np.linalg.norm(uList,axis=1))
    plt.plot([0,1],[0,1],'k--')
    plt.xlabel('Expected Std'); plt.ylabel('Abs(Error on VisCorr)'); plt.title('Noise Model Based on Beta and T_obs (color = u)')

    #%% Error vs. Data
    plt.figure(27); plt.clf()
    lincalScatter(np.abs(visCorrs), np.abs(uCal.computeErrors(betas, Sigmas, Ds)), color=np.linalg.norm(uList,axis=1), figNum=27)
    plt.xlabel('abs(visibility correlation)'); plt.ylabel('abs(error)'); plt.title('Bootstrapped Variance Model (color = u)')  
    plt.plot([0,1],[0,1],'k--')    
    
    plt.figure(28); plt.clf()
    lincalScatter(np.abs(visCorrs), np.abs(uCal.computeErrors(allResults[0][0], allResults[0][1], allResults[0][2])), color=np.linalg.norm(uList,axis=1), figNum=28)
    plt.xlabel('abs(visibility correlation)'); plt.ylabel('abs(error)'); plt.title('Noise Model Based on Beta and T_obs (color = u)')
    plt.plot([0,1],[0,1],'k--')
    
    #%% Error vs. Params
    betaIndex = {uCal.chans[i]: i for i in range(uCal.nChans)}
    SigmaIndex = {uCal.uBins[i]: i for i in range(uCal.nuBins)}
    DIndex = {uCal.duBins[i]: i for i in range(uCal.nduBins)}

    trial = 0
    figNum = 299
    
    beta1List = [allResults[trial][0][betaIndex[ch1]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
    beta2List = [allResults[trial][0][betaIndex[ch2]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
    betaProduct = [allResults[trial][0][betaIndex[ch1]]*np.conj(allResults[trial][0][betaIndex[ch2]]) for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
    SigmaList = [allResults[trial][1][SigmaIndex[entry['uBin']]] for entry in uCal.blChanPairs.values()]
    DList = [allResults[trial][2][DIndex[entry['duBin']]] for entry in uCal.blChanPairs.values()]
    samplesList = [entry['samples'] for entry in uCal.blChanPairs.values()]
    duList = np.asarray([entry['du'] for entry in uCal.blChanPairs.values()])
    uList = np.asarray([entry['u'] for entry in uCal.blChanPairs.values()])
    ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
    ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
    visCorrs = [entry['visCorr'] for entry in uCal.blChanPairs.values()]
    
    betaStds = np.std(np.abs(np.asarray(allResults[trial][4]).T),axis=1)
    SigmaStds = np.std(np.abs(np.asarray(allResults[trial][5]).T),axis=1)
    DStds = np.std(np.abs(np.asarray(allResults[trial][6]).T),axis=1)
    beta1StdList = [betaStds[betaIndex[ch1]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
    beta2StdList = [betaStds[betaIndex[ch2]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
    SigmaStdList = [SigmaStds[SigmaIndex[entry['uBin']]] for entry in uCal.blChanPairs.values()]
    DStdList = [DStds[DIndex[entry['duBin']]] for entry in uCal.blChanPairs.values()]
    absErrorList = np.abs(uCal.computeErrors(allResults[trial][0], allResults[trial][1], allResults[trial][2]))



    def errorComparisonPlot(x,y,color,figNum,xlabel,ylabel,title):
        plt.figure(figNum); plt.clf()
        lincalScatter(np.abs(x), np.abs(y), color=color, figNum=figNum)
        plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title('iteration ' + str(trial) + ': ' + title)
        plt.plot([0,1e6*np.median(x)],[0,1e6*np.median(y)],'k--')
        plt.colorbar()

#    figNum+=1; errorComparisonPlot(np.asarray(samplesList)**-.5, absErrorList, ch1List, figNum, '1/sqrt(samples)', 'abs(error)', 'color = chan1')

#    figNum+=1; errorComparisonPlot(beta1List, absErrorList, ch1List, figNum, 'beta1', 'abs(error)', 'color = chan1')
#    figNum+=1; errorComparisonPlot(beta2List, absErrorList, ch2List, figNum, 'beta2', 'abs(error)', 'color = chan2')
#    figNum+=1; errorComparisonPlot(betaProduct, absErrorList, np.abs(beta1List), figNum, 'beta1*beta2', 'abs(error)', 'color = beta1')
    figNum+=1; errorComparisonPlot(np.asarray(betaProduct) * np.asarray(samplesList)**-.5, absErrorList, np.abs(beta1List), figNum, 'beta1*beta2/sqrt(samples)', 'abs(error)', 'color = beta1')
    figNum+=1; errorComparisonPlot(SigmaList, absErrorList, np.linalg.norm(uList,axis=1), figNum, 'Sigma', 'abs(error)', 'color = u')
#    figNum+=1; errorComparisonPlot(DList, absErrorList, duList, figNum, 'Ds', 'abs(error)', 'color = du')

#    figNum+=1; errorComparisonPlot(ch1List, absErrorList, beta1List, figNum, 'chan1', 'abs(error)', 'color = beta1')
#    figNum+=1; errorComparisonPlot(ch2List, absErrorList, beta2List, figNum, 'chan2', 'abs(error)', 'color = beta2')
#    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), absErrorList, SigmaList, figNum, 'u', 'abs(error)', 'color = Sigma')
#    figNum+=1; errorComparisonPlot(duList, absErrorList, DList, figNum, 'du', 'abs(error)', 'color = D')
    
#    figNum+=1; errorComparisonPlot(beta1StdList, absErrorList, ch1List, figNum, 'beta1Std', 'abs(error)', 'color = chan1')
#    figNum+=1; errorComparisonPlot(beta1StdList, absErrorList, ch2List, figNum, 'beta1Std', 'abs(error)', 'color = chan2')
    figNum+=1; errorComparisonPlot(SigmaStdList, absErrorList, np.linalg.norm(uList,axis=1), figNum, 'SigmaStd', 'abs(error)', 'color = u')
#    figNum+=1; errorComparisonPlot(DStdList, absErrorList, duList, figNum, 'DStd', 'abs(error)', 'color = du')
    
#    figNum+=1; errorComparisonPlot(np.abs(np.asarray(ch1List)-np.asarray(ch2List)), absErrorList, np.linalg.norm(uList,axis=1), figNum, '|delta channel|', 'abs(error)', 'color = u')
    figNum+=1; errorComparisonPlot(np.abs(visCorrs), np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), np.linalg.norm(uList,axis=1), figNum, '|visCorr|', '|(error/visCorr)|', 'color = u')
    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), betaProduct, figNum, '|u|', '|(error/visCorr)|', 'color = betaProduct')
    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), np.asarray(betaProduct) * np.asarray(samplesList)**-.5, np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), figNum, '|u|', 'beta1*beta2/sqrt(samples)', 'color = |(error/visCorr)|')

    remainingError = np.asarray(absErrorList) / np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5)
    figNum+=1; errorComparisonPlot(SigmaList, remainingError, np.linalg.norm(uList,axis=1), figNum, 'SigmaList', 'abs(error)/(beta1*beta2*samples**-.5)', 'color = u')
    figNum+=1; errorComparisonPlot(SigmaStdList, remainingError, np.linalg.norm(uList,axis=1), figNum, 'SigmaStd', 'abs(error)/(beta1*beta2*samples**-.5)', 'color = u')
    


    #%%
    from scipy.stats import linregress
    toTest = {'1/sqrt(samples)': np.asarray(samplesList)**-.5,
            'beta1': np.abs(beta1List),
            'beta2': np.abs(beta2List),
            'betaProduct': np.abs(betaProduct),
            'betaProduct/samples**.5:': np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5),
            'deltaChan': np.abs(np.asarray(ch1List)-np.asarray(ch2List)),
            'Sigma': np.abs(SigmaList),
            'Ds': np.abs(DList),     
            'beta1Std': beta1StdList,
            'beta2Std': beta2StdList,
            'SigmaStd': SigmaStdList,
            'DStd': DStdList,
            'beta1Std*beta2Std': np.asarray(beta1StdList)*np.asarray(beta2StdList),
            'betaProduct*abs(Sigma)/samples**.5':  np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5) / np.asarray(np.abs(SigmaList)),
            'abs(visCorr)': np.abs(visCorrs),
            'u': np.linalg.norm(uList,axis=1),
            }

    #remainingError = absErrorList
#    remainingError = np.asarray(absErrorList) / np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5)
    remainingError = np.asarray(absErrorList)/np.asarray(np.abs(visCorrs))

    for testName, testVector in toTest.items():
        print '\n', testName
        slope, intercept, r_value, p_value, std_err = linregress(testVector,remainingError)
        print 'r =', r_value
        print 'slope =', slope
            
    
    
    
    
    
#%%
#     inferredBetas = {chan: [] for chan in uCal.chans}
#     betaDict = {uCal.chans[i]: betas[i] for i in range(uCal.nChans)}
#     for (ch1,bl1,ch2,bl2),entry in uCal.blChanPairs.items():
#         inferredBetas[ch1].append(entry['visCorr'] / betaDict[ch2].conj() / SigmaDict[entry['uBin']] / DDict[entry['duBin']])
#         inferredBetas[ch2].append(entry['visCorr'] / betaDict[ch1] / SigmaDict[entry['uBin']] / DDict[entry['duBin']])
#     inferredBetasStd = np.asarray([np.std(inferredBetas[chan])/(len(inferredBetas[chan])**.5) for chan in uCal.chans])
#        plt.plot(uCal.chans, inferredBetasStd)

#%%
#chanCompiledList = {chan: [] for chan in uCal.chans}
#chanCompiledList2 = {chan: [] for chan in uCal.chans}
#errorList = uCal.computeErrors(betas, Sigmas, Ds)
#visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
#modelList = visCorrs# - errorList
#ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
#ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
#    
#for f1,f2,error,Nii,model in zip(ch1List,ch2List,errorList,noiseCovDiag,modelList):
#    chanCompiledList[f1].append(error/model)#((2*Nii)**.5))
#    chanCompiledList[f2].append(error/model)#((2*Nii)**.5))
#    chanCompiledList2[f1].append(error/(2*Nii)**.5)
#    chanCompiledList2[f2].append(error/(2*Nii)**.5)
#chanAvgErrors = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList[chan]))) for chan in uCal.chans])
#chanAvgErrors2 = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList2[chan])))for chan in uCal.chans])
#plt.figure(4); plt.clf()
#plt.plot(uCal.chans, chanAvgErrors,'.')
#plt.plot(uCal.chans, chanAvgErrors2,'.-')
#plt.ylabel('Channel-Averaged, Noise Weighted Errors'); plt.xlabel('Channel')
##     badChans = np.asarray(uCal.chans)[chanAvgErrors > 2.5]
##     if len(badChans) > 0: print 'Channels with average sigma > 2.5: ', badChans
#
#print np.asarray(uCal.chans)[chanAvgErrors > .5]
#
#%%
#chanCompiledList = {chan: [] for chan in uCal.chans}
#for (ch1,bl1,ch2,bl2),entry in uCal.blChanPairs.items():
#    chanCompiledList[ch1].append(entry['samples'])
#    chanCompiledList[ch2].append(entry['samples'])
#chanTotalSamples = np.asarray([np.sum(np.asarray(chanCompiledList[chan])) for chan in uCal.chans])
#plt.figure(444); plt.clf()
#plt.plot(uCal.chans, chanTotalSamples,'.')
#print np.asarray(uCal.chans)[chanTotalSamples < np.mean(chanTotalSamples)/2]