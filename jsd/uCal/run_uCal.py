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
import uCal as uc

#############################################
#   Set-up script specific to PAPER 128. 
#   TODO: Hardcoded for now.
#############################################

regenerateEverything = False

if regenerateEverything:
    dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")]
    #dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.57058.xx.uvcRRE.npz")] 
    pol='xx'
    alwaysFlaggedChannels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 169, 183, 185, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 75, 76, 77]
    alsoFlagTheseChannels = [14, 15, 55, 101, 123, 179, 180, 181, 182, 184, 186]    
    flaggedChannels = sorted(list(set(alwaysFlaggedChannels).union(set(alsoFlagTheseChannels))))

    seps = np.arange(1,16) #* 15m  = baseline lengths
    dx = 15.0
    separations = dx*seps
    freqs = np.arange(.1,.2,.1/203)
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
        print 'Now generating baseline redundancy dictionary from calfile..'
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
    uReds = uc.uCalReds(freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau=.3) #just pass in freqs
    uReds.applyuCut(uMin=25, uMax=150)
    uReds.applyChannelFlagCut(flaggedChannels) 
    uCal = uc.uCalibrator(uReds.getBlChanPairs())
    uCal.computeVisibilityCorrelations(data, samples)
    pickle.dump([uReds, uCal], open('./Data/uCalData.p', 'wb'))
else:
    uReds, uCal = pickle.load(open('./Data/uCalData.p','rb'))

uCal.setupBinning(uBinSize = .72**.5, duBinSize = 5.0/203) #TODO: investigate uBinSize = .5 (nyquist sampling)
#uCal.setupBinning(uBinSize = .5, duBinSize = 5.0/203)

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
    noiseCovDiag = uCal.generateNoiseCovariance(betas) #updated each cycle based on improved result for beta

noiseCovDiag = uCal.renormalizeNoise(betas, Sigmas, Ds, noiseCovDiag)
betas, Sigmas, Ds, chiSqPerDoF = uCal.performLincalIteration(betas, Sigmas, Ds, noiseCovDiag, alpha = .5)    
print 'Final, noise-median-renormalized chi^2/dof = ' + str(chiSqPerDoF)

#TODO: add in polynomial fitting of final beta solution
#TODO: take out sky part from final solution

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
    ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
    ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])

    #%%Bandpass
    plt.figure(1); plt.clf()
    inferredErrorsOnAbsBeta = np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    plt.errorbar(np.arange(.1,.2,.1/203)[uCal.chans],np.abs(betas),yerr=inferredErrorsOnAbsBeta)

    plt.xlabel('Frequency (GHz)'); plt.ylabel('Abs(Lincal Bandpass)');

    #%%Predicted vs. Observed Scatter
    plt.figure(2); plt.clf()
    visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
    predictedCorrs = visCorrs - uCal.computeErrors(betas, Sigmas, Ds)
    lincalScatter(np.abs(predictedCorrs), np.abs(visCorrs), color=duList, figNum=2, ys='log', xs='log', title = '')
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
        chanCompiledList[f1].append(error/(2*Nii**.5))
        chanCompiledList[f2].append(error/(2*Nii**.5))
    chanAvgErrors = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList[chan]))) for chan in uCal.chans])
    plt.figure(4); plt.clf()
    plt.plot(uCal.chans, chanAvgErrors,'.')
    plt.ylabel('Channel-Averaged, Noise Weighted Errors'); plt.xlabel('Channel')
    badChans = np.asarray(uCal.chans)[chanAvgErrors > 3]
    if len(badChans) > 0: print 'Channels with average sigma > 3: ', badChans

    plt.show()