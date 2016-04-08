import numpy as np
from scipy.sparse import csr_matrix
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import scipy
import omnical
import time
import aipy as a
import optparse, sys, os
import capo

#############################################
#   uCal Core Classes
#############################################

#TODO needs testing on a 2D or 3D array. In particular, there's a question of whether the code that skips entries to make the loop fasterever misses a good pair
class uCalReds():
    """This class takes a list frequencies (in GHz), a list of baselines, anda  way to convert them into numpy arrays in meters.
    It saves a dictionary that maps (ch1,bl1,ch2,bl2) to u and deltau. 
    Physical separations must all have the same orentiation, otherwise this won't work.
    Only includes baselines-frequency pairs obeying the Deltau threshold (default 0.3)."""

    #TODO: remove chan2FreqDict
    def __init__(self, freqs, bls, chan2FreqDict, bl2SepDict, maxDeltau = .3, verbose=False):
        self.verbose = verbose
        if self.verbose: print 'Now finding all baseline/frequency pairs...'
        self.maxDeltau = maxDeltau
        chans = range(len(freqs))
        self.blChanPairs = {}
        freqChanPairs = sorted([[chan2FreqDict[chan], chan] for chan in chans])
        sepBLPairs = sorted([[np.asarray(bl2SepDict[bl]), bl] for bl in bls], key=lambda pair: np.linalg.norm(pair[0])) #sorts by length
        for i,(f1,ch1) in enumerate(freqChanPairs):
            for f2,ch2 in freqChanPairs[i+1:]:
                for j,(sep2,bl2) in enumerate(sepBLPairs):
                    for sep1,bl1 in sepBLPairs[j+1:]:
                        deltau = np.linalg.norm(f1*sep1 - f2*sep2) 
                        if deltau < maxDeltau and not self.blChanPairs.has_key((ch1,bl1,ch2,bl2)) and not self.blChanPairs.has_key((ch2,bl2,ch1,bl1)):
                            u = (f1*sep1 + f2*sep2)/2.0 
                            self.blChanPairs[(ch1,bl1,ch2,bl2)] = (u, deltau)
        if self.verbose: print "    " + str(len(self.blChanPairs)) + " baseline/frequency pairs identified with delta u < " + str(maxDeltau)

    def applyuCut(self, uMin=25, uMax=150):
        for key,value in self.blChanPairs.items():
            if np.linalg.norm(value[0]) < uMin or np.linalg.norm(value[0]) > uMax: del[self.blChanPairs[key]]
        if self.verbose: print "    " + str(len(self.blChanPairs)) + " baseline/frequency pairs remain after requiring " + str(uMin) + " < u < " + str(uMax)

    def applyChannelFlagCut(self, flaggedChannels):
        for (ch1,bl1,ch2,bl2) in self.blChanPairs.keys():
            if ch1 in flaggedChannels or ch2 in flaggedChannels: del[self.blChanPairs[(ch1,bl1,ch2,bl2)]]
        if self.verbose: print "    " + str(len(self.blChanPairs)) + " baseline/frequency pairs remain after flagging " + str(len(set(flaggedChannels))) + " channels."

    def getBlChanPairs(self): return self.blChanPairs


class uCalibrator():
    """This class contains the main routines for performing uCal. 
    It is intialized with a dictionary of baseline-channels pairs (see uCalReds class for details) 
    and optionally u and du bin sizes. The constructor figures out the relevant binning."""

    def __init__(self, blChanPairs):
        self.lincalIterations = 0
        #Internal format for blChanPairs that will be populated with data and binning info
        self.blChanPairs = {key: {'u': u, 'du': du} for key,(u,du) in blChanPairs.items()}
        self.nPairs = len(self.blChanPairs)
        self.binningIsSetup = False
        self.visibilitiesAreCorrelated = False

    def getBlChanPairs(self): return self.blChanPairs

    def computeVisibilityCorrelations(self, data, samples):
        """This function computes visibility correlations from data dictionaries and samples dictionaries (which reflects flags and redundancies).
        These dictionaries must be in the standard PAPER format. The results are stored inside self.blChanPairs with keys 'visCorr' and 'samples'."""
        for (ch1,bl1,ch2,bl2) in self.blChanPairs.keys():
            w = np.logical_and(samples[bl1][:,ch1] != 0, samples[bl2][:,ch2] != 0)
            if np.all(np.logical_not(w)):
                self.blChanPairs[(ch1,bl1,ch2,bl2)]['visCorr'] = 0.0+0.0j
                self.blChanPairs[(ch1,bl1,ch2,bl2)]['samples'] = 0
            else: 
                self.blChanPairs[(ch1,bl1,ch2,bl2)]['visCorr'] = np.average((data[bl1][:,ch1]*np.conj(data[bl2][:,ch2]))[w])
                self.blChanPairs[(ch1,bl1,ch2,bl2)]['samples'] = np.sum((samples[bl1][:,ch1] * samples[bl2][:,ch2])[w])
        self.visibilitiesAreCorrelated = True

    #TODO: need to figure out what to do if there are additional flagged channels. I guess, in general, I need to figure out what to do about flagged channels.


    def setupBinning(self, uBinSize = .72**-.5, duBinSize = 5.0/203):
        """Given a size of each bin in u (in wavelengths) and each bin in Delta u (in wavelengths), this function initializes the proper """
        #Find unique chans and count chans.
        self.chans = sorted(np.unique([[ch1,ch2] for (ch1,bl1,ch2,bl2) in self.blChanPairs.keys()]))
        self.nChans = len(self.chans)
        #Determine binning: first assign integers to u and du
        for key,value in self.blChanPairs.items():
            self.blChanPairs[key]['uBin'] = tuple(np.floor(1.0*(value['u'])/uBinSize).astype(int).tolist())
            self.blChanPairs[key]['duBin'] = np.floor(1.0*(value['du'])/duBinSize).astype(int)
        #Now find and sort the unique values of those integers
        self.uBins = sorted({value['uBin']: None for key,value in self.blChanPairs.items()}.keys(), key=lambda uBin: np.linalg.norm(np.asarray(uBin)))
        self.duBins = sorted({value['duBin']: None for key,value in self.blChanPairs.items()}.keys())
        #Now determine bin centers
        self.uBinCenters = {uBin: [] for uBin in self.uBins}
        self.duBinCenters = {duBin: [] for duBin in self.duBins}
        for key,value in self.blChanPairs.items():
            self.uBinCenters[value['uBin']].append(value['u'])
            self.duBinCenters[value['duBin']].append(value['du'])
        self.uBinCenters = [np.mean(np.asarray(self.uBinCenters[uBin]), axis = 0) for uBin in self.uBins]
        self.nuBins = len(self.uBinCenters)
        self.duBinCenters = [np.mean(self.duBinCenters[duBin]) for duBin in self.duBins]
        self.nduBins = len(self.duBinCenters)
        self.binningIsSetup = True

    def normalize(self, betas, Sigmas, Ds):
        """This function normalizes betas, Sigmas, and Ds according to the following convention:\n
            - Mean absolute value of betas = 1.0
            - Sum of phases of beta = 0.0
            - Max value of Ds = 1.0 """
        bandpassMean = np.average(np.abs(betas))
        bandpassAngleMean = np.average(np.angle(betas))
        DsMax = np.max(np.abs(Ds))
        return betas / bandpassMean * np.exp(-1.0j * bandpassAngleMean), Sigmas * bandpassMean**2 * DsMax, Ds / DsMax

    def performLogcal(self):
        """This function returns the logcal result beta, Sigma, and D in the order specified."""
        if not self.binningIsSetup: raise RuntimeError('Binn')
        if not self.visibilitiesAreCorrelated: raise RuntimeError('')
        #Construct logcal matrices
        chan2Col = {chan: col for col,chan in zip(range(0,self.nChans),self.chans)}
        uBin2Col = {uBin: col for col,uBin in zip(range(self.nChans,self.nChans+self.nuBins),self.uBins)}
        duBin2Col = {duBin: col for col,duBin in zip(range(self.nChans+self.nuBins,self.nChans+self.nuBins+self.nduBins),self.duBins)}
        Acoeffs, Bcoeffs, rowIndices, colIndices = np.zeros(self.nPairs*4), np.zeros(self.nPairs*4), np.zeros(self.nPairs*4), np.zeros(self.nPairs*4)
        for n,((ch1,bl1,ch2,bl2),entry) in enumerate(self.blChanPairs.items()):
            rowIndices[4*n:4*n+4] = n
            colIndices[4*n:4*n+4] = [chan2Col[ch1], chan2Col[ch2], uBin2Col[entry['uBin']], duBin2Col[entry['duBin']]]
            Acoeffs[4*n:4*n+4] = [1.0, 1.0, 1.0, 1.0]
            Bcoeffs[4*n:4*n+4] = [1.0, -1.0, 1.0, 1.0]

        self.logcalA, self.logcalB = [csr_matrix((coeffs,(rowIndices,colIndices))) for coeffs in (Acoeffs,Bcoeffs)]
        self.logcalAtA, self.logcalBtB = [M.conjugate().transpose().dot(M).toarray() for M in (self.logcalA, self.logcalB)]

        logcalZeroEVs = [len(MtM) - np.linalg.matrix_rank(MtM) for MtM in (self.logcalAtA, self.logcalBtB)]
        if not logcalZeroEVs == [2,2]: print "    WARNING: Logcal's AtA (real part) has " + str(logcalZeroEVs[0]) + " zero eigenvalues. Logcal's BtB (imag part) has " + str(logcalZeroEVs[0]) + " zero eigenvalues. They should both have 2 each."

        #Perform logcal
        y = np.asarray([np.log(value['visCorr']) for value in self.blChanPairs.values()])
        xhatReal = np.linalg.pinv(self.logcalAtA).dot((self.logcalA.conjugate().T).dot(np.real(y)))
        xhatImag = np.linalg.pinv(self.logcalBtB).dot((self.logcalB.conjugate().T).dot(np.imag(y)))
        result = np.exp(xhatReal + 1.0j*xhatImag)
        betas, Sigmas, Ds = self.normalize(result[0:self.nChans], result[self.nChans: self.nChans+self.nuBins], result[self.nChans+self.nuBins:self.nChans+self.nuBins+self.nduBins])
        return betas, Sigmas, Ds


    def generateNoiseCovariance(self, betas):
        """Creates a noise covariance that is proportional to beta(ch1)**2 * beta(ch2)**2 / nSamples where
        nSamples is nBaselines_1 * nBaselines_2 * nIntegrations. Overall scaling is arbitrary."""
        betaAbsDict = {self.chans[n]: np.abs(betas[n]) for n in range(self.nChans)}
        noiseCovDiag = np.asarray([(betaAbsDict[ch1]**2 * betaAbsDict[ch2]**2) / (1.0 * entry['samples']) for (ch1,bl1,ch2,bl2),entry in self.blChanPairs.items()])
        noiseCovDiag = noiseCovDiag / np.median(noiseCovDiag) * 4e-10 #approximate renormalization
        return noiseCovDiag

    def uCalDicts(self, betas, Sigmas, Ds):
        """This function turns beta, Sigma, and D into dictionaries indexed by their bin/channel number."""
        betaDict = {self.chans[n]: betas[n] for n in range(self.nChans)}
        SigmaDict = {self.uBins[n]: Sigmas[n] for n in range(self.nuBins)}
        DDict = {self.duBins[n]: Ds[n] for n in range(self.nduBins)}
        return betaDict, SigmaDict, DDict

    def computeErrors(self, betas, Sigmas, Ds):
        """This function computes the difference between the measurements and the predicted visibility correlations (in the standard order of self.blChanPairs.items())."""
        betaDict, SigmaDict, DDict = self.uCalDicts(betas, Sigmas, Ds)
        return np.asarray([entry['visCorr'] - betaDict[ch1]*np.conj(betaDict[ch2])*SigmaDict[entry['uBin']]*DDict[entry['duBin']] for (ch1,bl1,ch2,bl2),entry in self.blChanPairs.items()])

    def performLincalIteration(self, betas, Sigmas, Ds, noiseCovDiag, alpha = .5):
        """This function performs lincal. It takes:\n
            - A starting guess for beta (use logcal to initialize this the first them, then perform iteratively)
            - A starting guess for Sigma (to be improved iteratively)
            - A starting guess for D (to be improved iteratively)
            - A model for the noise covariance diagonal (in the standard order of self.blChanPairs.items())
            - Optionally, alpha, an amount by which to update guesses, between 0 and 1.
                -Lower alpha converges slower.
                -Larger alpha is more likely to diverge."""
        self.lincalIterations += 1
        chan2Col = {chan: col for col,chan in zip(range(0, 2*self.nChans, 2), self.chans)}
        uBin2Col = {uBin: col for col,uBin in zip(range(2*self.nChans, 2*self.nChans+2*self.nuBins, 2), self.uBins)}
        duBin2Col = {duBin: col for col,duBin in zip(range(2*self.nChans+2*self.nuBins, 2*self.nChans+2*self.nuBins+2*self.nduBins, 2), self.duBins)}
        betaDict, SigmaDict, DDict = self.uCalDicts(betas, Sigmas, Ds)

        coeffs, rowIndices, colIndices = np.zeros(self.nPairs*16), np.zeros(self.nPairs*16), np.zeros(self.nPairs*16)
        for n,((ch1,bl1,ch2,bl2),entry) in enumerate(self.blChanPairs.items()):
            bbstarSD = betaDict[ch1]*np.conj(betaDict[ch2])*SigmaDict[entry['uBin']]*DDict[entry['duBin']] #(beta)(beta^*)(Sigma)(D)
            bbstarS = betaDict[ch1]*np.conj(betaDict[ch2])*SigmaDict[entry['uBin']] #(beta)(beta^*)(Sigma)
            bbstarD = betaDict[ch1]*np.conj(betaDict[ch2])*DDict[entry['duBin']] #(beta)(beta^*)(D)
            ch1Col, ch2Col, uCol, duCol = chan2Col[ch1], chan2Col[ch2], uBin2Col[entry['uBin']], duBin2Col[entry['duBin']]
            rowIndices[16*n:18*n+6] = n #the first 8 terms are on this row
            rowIndices[16*n+8:16*n+16] = n+self.nPairs #the next 8 terms are on the corresponding imaginary part row
            for i,colIndex in enumerate([ch1Col,ch2Col,ch1Col+1,ch2Col+1,uCol,uCol+1,duCol,duCol+1,ch1Col,ch2Col,ch1Col+1,ch2Col+1,uCol,uCol+1,duCol,duCol+1]):
                colIndices[16*n+i] = colIndex #these are eta1, eta2, phi1, phi2, sigma, psi, etc.
            coeffList = [np.real(bbstarSD), np.real(bbstarSD), -np.imag(bbstarSD), np.imag(bbstarSD), np.real(bbstarD), -np.imag(bbstarD), np.real(bbstarS), -np.imag(bbstarS),
                         np.imag(bbstarSD), np.imag(bbstarSD), np.real(bbstarSD), -np.real(bbstarSD), np.imag(bbstarD), np.real(bbstarD), np.imag(bbstarS), np.real(bbstarS)]
            for i,coeff in enumerate(coeffList): 
                coeffs[16*n+i] = coeff #these are the coefficients of those terms
        
        self.A = csr_matrix((coeffs,(rowIndices,colIndices)))
        Ninv = csr_matrix((np.append((noiseCovDiag)**-1,(noiseCovDiag)**-1), (np.arange(2*self.nPairs), np.arange(2*self.nPairs))))
        self.AtNinvA = (self.A.conjugate().transpose().dot(Ninv)).dot(self.A).toarray()
        if self.lincalIterations == 1: 
            lincalZeroEVs = len(self.AtNinvA) - np.linalg.matrix_rank(self.AtNinvA)
            if not lincalZeroEVs == 4: print "    WARNING: Lincal's AtNinvA has " + str(lincalZeroEVs) + " zero eigenvalues. It should have 4."

        deltas = self.computeErrors(betas, Sigmas, Ds)
        xHat = np.linalg.pinv(self.AtNinvA).dot(self.A.T.conjugate().dot(Ninv.dot(np.append(np.real(deltas),np.imag(deltas)))))
        newBetas = np.asarray([betaDict[chan]*(1+alpha*(xHat[chan2Col[chan]] + 1.0j*xHat[chan2Col[chan]+1])) for chan in self.chans])
        newSigmas = np.asarray([SigmaDict[uBin] + alpha*(xHat[uBin2Col[uBin]] + 1.0j*xHat[uBin2Col[uBin]+1]) for uBin in self.uBins])
        newDs = np.asarray([DDict[duBin] + alpha*(xHat[duBin2Col[duBin]] + 1.0j*xHat[duBin2Col[duBin]+1]) for duBin in self.duBins])
        betas, Sigmas, Ds = self.normalize(newBetas, newSigmas, newDs)

        errors = self.computeErrors(betas, Sigmas, Ds)
        chiSqPerDoF = np.average(np.append(np.real(errors)**2, np.imag(errors)**2) * Ninv)
        return betas, Sigmas, Ds, chiSqPerDoF

    def renormalizeNoise(self, betas, Sigmas, Ds, noiseCovDiag):
        """Returns a new noise covariance diagonal that has been renormalized such that the median observed error reflects the median noise covariance. """
        errors = self.computeErrors(betas, Sigmas, Ds)
        return noiseCovDiag * (np.median(np.real(errors)**2)+np.median(np.imag(errors)**2)) / (2*np.median(noiseCovDiag))

def save2npz(filename, dataFiles, bandpass, weights, betas, Sigmas, Ds, noiseCovDiag, A):
    """Saves necessary information to use this uCal solution on data or to combine with with other times:\n
        - filename: output filename
        - dataFiles: list of files that went into this uCal solution
        - bandpass: best guess of the bandpass (sky power law corrected), length = len(freqs)
        - weights: 0 for flagged channels, 1 otherwise, length = len(freqs)
        - betas, Sigmas, Ds: raw result of uCal
        - noiseCovDiag: noise on the real (or imaginary) part of the visibiltiy correlations
        - A: lincal A matrix. """
    result = {}
    result['dataFiles'] = dataFiles
    result['bandpass'] = bandpass
    result['weights'] = weights
    result['betas'] = betas
    result['Sigmas'] = Sigmas
    result['Ds'] = Ds
    result['A'] = A
    result['noiseCovDiag'] = noiseCovDiag
    np.savez(filename, **result)

def loadAndCombine(fileList):
    invCovList = []
    invCovWeightedResultList = []
    for file in fileList:
        result = np.load(file)
        

# if __name__ == '__main__':
#     loadAndCombine(['/Users/jsdillon/Desktop/capo/jsd/uCal/Data/zen.2456943.65409.xx.uvcRRE.uCal.npz'])






