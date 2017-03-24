import numpy as np
#from autograd import numpy as np
#import autograd
#from autograd.convenience_wrappers import hessian_vector_product as hvp
from scipy.sparse import csr_matrix
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import scipy
import omnical
import time
import aipy as a
import optparse, sys, os
import capo
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from itertools import chain
import copy

#############################################
#   uCal Core Class
#############################################

class uCalibrator():
    """This is the core class for performing uCal. It stores the relevant visibilities and information about them. It also performs logcal and lincal."""

#     ###########################################################################
#     #   Trying to use python's built in stuff
#     ###########################################################################


#     def objectiveFunc(self, allParams):
#         allParams = np.array(allParams[0:len(allParams)/2]) + 1.0j*np.array(allParams[len(allParams)/2:])
#         betas = allParams[0:self.nChans]
#         Sigmas = allParams[self.nChans:]
#         sumSquareErrors = 0
#         for n,entry in enumerate(self.blChanPairs.values()):
#             sumSquareErrors = sumSquareErrors + np.abs(entry['vis'] - self.bandpassFunc(entry, betas) * self.uIntegralFunc(entry, Sigmas))**2
        
#         return sumSquareErrors


#     def renormalize2(self, allParams, absMeanGoal = 1.0):
#         """Sets the mean phase of the betas to 0 and the mean amplitude to 1."""
#         betas = allParams[0:self.nChans]
#         Sigmas = allParams[self.nChans:]

#         betaAbsMean = np.mean(np.abs(betas))
#         betaAngleMean = np.mean(np.angle(betas))
#         newBetas = np.array(betas) * absMeanGoal/betaAbsMean * np.exp(-1.0j * betaAngleMean)
#         newSigmas = np.array(Sigmas)
#         newSigmas[0:self.nuSamples*(self.skyFreqOrder+1)] = newSigmas[0:self.nuSamples*(self.skyFreqOrder+1)] * betaAbsMean/absMeanGoal / np.exp(-1.0j * betaAngleMean)
#         #newSigmas[self.beamWidthIndices[0]:self.beamWidthIndices[-1]+1] = np.real(np.array(newSigmas[self.beamWidthIndices[0]:self.beamWidthIndices[-1]+1]))
#         #TODO: figure out whether there's a better way to handle this.
#         return np.concatenate((newBetas, newSigmas))

#     def callbackFunc(self, allParams):
#         print np.mean(np.abs(allParams[0:self.nChans]))
#         #allParams = self.renormalize2(allParams)
#         print self.objectiveFunc(allParams)

#     def minimizeError(self, betas, Sigmas):
#         print 'Now minimizing error...'
#         from scipy.optimize import minimize

#         g_auto = autograd.grad(self.objectiveFunc)
#         hessvp = hvp(self.objectiveFunc)
#         #cons = ({'type': 'eq', 'fun': lambda x: np.mean(np.abs(x[0:self.nChans]))-1.0},
# #                        {'type': 'eq', 'fun': lambda x: np.mean(np.angle(x[0:self.nChans]))})
#         x0 = np.concatenate((np.real(betas), np.real(Sigmas), np.imag(betas), np.imag(Sigmas)))
#         result = minimize(self.objectiveFunc, x0, method='trust-ncg', jac=g_auto, hessp =hessvp, callback=self.callbackFunc)#, constraints=cons)
#         #result = minimize(self.objectiveFunc, x0, method='BFGS', jac=g_auto, callback=self.callbackFunc)#, constraints=cons)
#         print result
#         return np.array(result[0][0:len(result[0])/2]) + 1.0j*np.array(result[0][len(result[0])/2:])



    ###########################################################################
    #   Specialized Functions that Depend on Model Paramterization
    ###########################################################################

    def bandpassFunc(self, blChanPair, betas, derivative=False):
        """This returns elements of the bandpass function beta, or derivaties and indices."""
        if derivative: return np.array([1.0]), [self.chanIndices[blChanPair['chan']]]
        else: return betas[self.chanIndices[blChanPair['chan']]]

    # TODO LIST:
    # 0) Write rudimentary post-logcal script DONE
    # 1) Write pre-computation script for each lincal step DONE
    # 2) Write uIntegralFunction for derivative=False DONE
    # 3) Update post-logcal script with beam fitting DONE
    # 4) Write uIntegralFunction with derivative = True

    def computeGaussianBeam(self, freq, fullRangeDeltaus):
        return np.exp(-(np.array(fullRangeDeltaus))**2 / (2*self.beamWidthDict[freq]**2))

    def FourierInterpolate(self, toInterp, pad):
        """Takes some evenly spaced function and returns a Fourier interpolation of that function at the padding level."""
        toPad = np.fft.fftshift(np.fft.fft(toInterp))
        padded = np.concatenate((np.zeros(int(np.ceil(len(toInterp)*pad+.5))), toPad, np.zeros(int(np.floor(len(toInterp)*pad-.5)))))
        interpolated = np.roll((2*pad+1)*np.fft.ifft(np.fft.fftshift(padded)),int(np.floor(pad)))
        return interpolated

    def FourierDownsample(self, toDownSample, pad):
        """Performs the inverse operation of FourierInterpolate. Used for calculating derivaties of uIntegralFunc."""
        finalLength = int(np.round(len(toDownSample) / (2*pad + 1)))
        shifted = np.roll(toDownSample, -int(np.floor(pad))) / (2*pad+1)
        FTed = np.fft.ifftshift(np.fft.fft(shifted))
        cropped = FTed[int(np.ceil(finalLength*pad+.5)): int(np.ceil(finalLength*pad+.5))+finalLength]
        return np.fft.ifft(np.fft.ifftshift(cropped)) 

    def uIntegralFunc(self, blChanPair, uIntegralParams, derivative=False):
        """This function returns a beam-weighted integral over a portion of the uv-plane (line, really). 
        It can also return derivatives with respect to the paramters that describe that plane or the beam."""

        #TODO: figure out delta u as the length element of the integral 

        u,chan = blChanPair['u'], blChanPair['chan']
        coarseuBin, fineuBin = blChanPair['coarseuBin'], blChanPair['fineuBin']
        beamCoarseIndices = blChanPair['beamCoarseIndices']

        freq = self.freqs[chan]
        coarseu = self.coarseus[coarseuBin]
        fineOffset = int(np.round((u - coarseu)/self.fineDu))
        beamWidth = self.beamWidthDict[freq]
        nFine = len(self.fineDeltaus)

        coarseSky = np.zeros(self.beamCoarseBins,dtype=complex)
        for order in range(self.skyFreqOrder+1): coarseSky += uIntegralParams[blChanPair['coarseSkyIndices'][order,:]] * blChanPair['coarseSkyBasisFuncs'][order,:]
        interpolatedSky = self.FourierInterpolate(coarseSky, self.padding)
        gaussBeam = self.gaussBeamDict[freq][self.beamFineBins/2-nFine/2 - fineOffset: self.beamFineBins/2 + nFine/2 - fineOffset]
        normalization = np.linalg.norm(gaussBeam) * self.fineDu

        #This doesn't seem like it's working. Maybe we need see if the analytic derivatives are matching numerical ones. 
        #ALso, I should check that nothing else is in a "testing" mode that is screwing with the results. 

        if derivative: 
            downsampledBeam = self.FourierDownsample(gaussBeam, self.padding)
            sampleDerivs = blChanPair['uParamBasisFuncs'] * downsampledBeam[blChanPair['uParamBeamIndices']] * normalization

            uSecondMoment = np.dot(interpolatedSky * gaussBeam, blChanPair['fineDeltausHere']**2) 
            beamWidthDerivs = np.array([(freq/.145)**n * normalization / beamWidth**3 * uSecondMoment for n in range(self.beamWidthOrder+1)])

            indices = blChanPair['uParamIndicesHere'] + self.beamWidthIndices
            derivs = np.concatenate((sampleDerivs, beamWidthDerivs))

            return derivs, indices
        return normalization * np.dot(gaussBeam, interpolatedSky)





        # self.interpolatedSky = interpolatedSky
        # self.coarseSky = coarseSky
        # 
        # print 'u:', u
        # print 'nearest coarse u:', coarseu
        # print 'fine center offset:', fineOffset
        # print 'full length of gaussBeam:', len(self.gaussBeamDict[freq])
        # print 'gauss beam extraction indices:', self.beamFineBins/2-nFine/2 - fineOffset, self.beamFineBins/2 + nFine/2 - fineOffset
        # print 'length of extracted gaussBeam', len(gaussBeam)
        # print 'len(gaussBeam) vs. len(interpolatedSky):', len(gaussBeam), len(interpolatedSky)
        # print 'index of fineDeltau = 0:', np.argwhere(np.abs(fineDeltausHere) < 1e-10)
        # print 'gaussBeam at fineDeltau = 0:', gaussBeam[np.argwhere(np.abs(fineDeltausHere) < 1e-10)[0]]
        # print 'this should give 0:', (self.fineDeltaus + coarseu)[np.argwhere(np.abs(fineDeltausHere) < 1e-10)[0]] - u



    def uIntegralPrecomputation(self, uIntegralParams):
        """This function is meant to be run before each iteration of lincal and does useful pre-computation for uIntegralFunc()"""
        self.beamWidthDict = {freq: np.dot(uIntegralParams[self.beamWidthIndices], np.array([(freq/.145)**n for n in range(self.beamWidthOrder+1)])) for freq in self.freqs}
        self.gaussBeamDict = {freq: self.computeGaussianBeam(freq, self.fineDu*np.arange(-self.beamFineBins/2.0+1, self.beamFineBins/2.0+1)) for freq in self.freqs}


    def renormalize(self, factorFunctions, factorParams, absMeanGoal = 1.0):
        """Sets the mean phase of the betas to 0 and the mean amplitude to 1."""
        betaAbsMean = np.mean(np.abs(factorParams[0]))
        betaAngleMean = np.mean(np.angle(factorParams[0]))
        newBetas = np.array(factorParams[0]) * absMeanGoal/betaAbsMean * np.exp(-1.0j * betaAngleMean)
        newSigmas = np.array(factorParams[1])
        newSigmas[self.uParamIndexList] *= betaAbsMean/absMeanGoal / np.exp(-1.0j * betaAngleMean)
        return [newBetas, newSigmas]

    def determineuParameterization(self):
        """This function figures out the maximum order in basis functions that any given u mode can support.
        It also creates a dictionary to map u,order -> index in Sigma."""
        blRange = {} #Figure out the largest and smallest us for each baseline
        for (bl,chan),entry in self.blChanPairs.items():
            if bl in blRange: blRange[bl] = [min(blRange[bl][0], entry['u']), max(blRange[bl][1], entry['u'])]
            else: blRange[bl] = [entry['u'], entry['u']]
        
        #Count the number of baselines that overlap a given u sample, then set self.uOrders based on that
        allbls = list(set([bl for (bl,chan) in self.blChanPairs.keys()]))
        self.uOrders = {}
        for u in self.coarseus:
            nOverlappingBaselines = np.sum([(u >= blRange[bl][0])*(u <= blRange[bl][1]) for bl in allbls])
            self.uOrders[u] = min(max(0, nOverlappingBaselines), self.skyFreqOrder)

        #Build up a dictionary that maps u and order to uParamIndex
        currentInd = 0
        self.uParamIndices = {}
        for order in range(self.skyFreqOrder+1):
            for u in self.coarseus:
                if self.uOrders[u] >= order:
                    self.uParamIndices[(u,order)] = currentInd
                    currentInd += 1;
        if self.verbose: print 'For lincal, we are trying to fit', currentInd, 'u sample parameters.'
        self.uParamIndexList = range(currentInd)
        self.beamWidthIndices = range(self.uParamIndexList[-1]+1, self.uParamIndexList[-1]+1 + self.beamWidthOrder + 1)

        #Preform some precomputation that can be used for every iteration 
        for entry in self.blChanPairs.values():
            u, chan, coarseuBin, fineuBin = entry['u'], entry['chan'], entry['coarseuBin'], entry['fineuBin']
            fineOffset = int(np.round((u - self.coarseus[coarseuBin])/self.fineDu))
            entry['beamCoarseIndices'] = np.arange(coarseuBin - (self.beamCoarseBins-1)/2,coarseuBin + (self.beamCoarseBins+1)/2)
            entry['uParamIndicesHere'] = [self.uParamIndices[(uHere,order)] for uHere in self.coarseus[entry['beamCoarseIndices']] for order in range(self.uOrders[uHere]+1)]
            entry['uParamOrders'] = [order for uHere in self.coarseus[entry['beamCoarseIndices']] for order in range(self.uOrders[uHere]+1)]
            entry['uParamBasisFuncs'] = [self.skyFreqBasis[chan,order] for uHere in self.coarseus[entry['beamCoarseIndices']] for order in range(self.uOrders[uHere]+1)]
            entry['uParamBeamIndices'] = [index for index, uHere in enumerate(self.coarseus[entry['beamCoarseIndices']]) for order in range(self.uOrders[uHere]+1)]
            entry['fineDeltausHere']  = self.fineDeltaus - fineOffset * self.fineDu 

            entry['coarseSkyIndices'] = np.zeros((self.skyFreqOrder+1, self.beamCoarseBins),dtype=int)
            entry['coarseSkyBasisFuncs'] = np.zeros((self.skyFreqOrder+1, self.beamCoarseBins))
            for uBin, uHere in enumerate(self.coarseus[entry['beamCoarseIndices']]):
                for order in range(self.skyFreqOrder+1):
                    if order < self.uOrders[uHere]+1:
                        entry['coarseSkyIndices'][order,uBin] = self.uParamIndices[(uHere,order)]
                        entry['coarseSkyBasisFuncs'][order,uBin] = self.skyFreqBasis[chan,order]

    def errfunc(self, beamParams, betas, uPlaneParams):
        uIntegralParams = np.concatenate((uPlaneParams,beamParams))
        self.uIntegralPrecomputation(uIntegralParams)
        errors = np.zeros(self.nVis, dtype=complex)
        for n,entry in enumerate(self.blChanPairs.values()):
            errors[n] = entry['vis'] - self.bandpassFunc(entry, betas) * self.uIntegralFunc(entry, uIntegralParams)
        return np.array(np.abs(errors), dtype='float')

        
    def initializeLogcalResults(self, allParams):
        self.determineuParameterization()

        betas = allParams[0:self.nChans]
        Sigmas = np.zeros(len(self.uParamIndexList), dtype=complex)
        Sigmas[0:self.nCoarseBins] = allParams[self.nChans:] / np.mean(self.skyFreqBasis[:,0]) / self.coarseDu

        if self.verbose: print 'Now initializing the frequency-dependent primary beam kernel by fitting...'
        intialBeamParamsGuess = [self.beamFWHMguess / (8*np.log(2))**.5] + [0.0] * self.beamWidthOrder
        bestFitBeamParams = np.asarray(leastsq(self.errfunc, intialBeamParamsGuess, args=(betas,Sigmas), xtol=.01)[0], dtype=complex)
        Sigmas = np.append(Sigmas, bestFitBeamParams)
        self.beamWidths = [np.dot(bestFitBeamParams, [(self.freqs[chan]/(.145+0.0j))**n for n in range(self.beamWidthOrder+1)]) for chan in self.chans]

        
        # #If the guessed beam width doesn't work (produces negative beam sizes), try deltau. If that fails, try deltau*2, deltau*3, etc.
        # for beamMult in range(1, int(self.beamFullRange / self.deltau)+1):
        #     intialBeamParamsGuess = [self.beamFWHMguess / (8*np.log(2))**.5] + [0.0] * self.beamWidthOrder
        #     bestFitBeamParams = np.asarray(leastsq(self.errfunc, intialBeamParamsGuess, args=(betas,Sigmas), xtol=.01)[0], dtype=complex)
        #     beamWidths = [np.dot(bestFitBeamParams, [(self.freqs[chan]/(.145+0.0j))**n for n in range(self.beamWidthOrder+1)]) for chan in self.chans]
        #     if np.min(beamWidths) > 0 or beamMult == int(self.beamFullRange / self.deltau): 
        #         Sigmas = np.append(Sigmas, bestFitBeamParams)
        #         break
        #     else: 
        #         print '    Unable to initialize beam with guessed width ' + str(self.beamFWHMguess) + ', trying ' + str(self.deltau*beamMult) + ' instead.'
        #         self.beamFWHMguess = self.deltau*beamMult

        return self.renormalize(None, [betas, Sigmas])

    def generateNoiseVariance(self, betas): 
        """Creates a model for the noise variance, which is assumed proportional to (beta^2) / (baseline redundancy)."""
        noiseVariance = np.asarray([np.abs(betas[self.chanIndices[chan]])**2 / self.redundancyDict[bl] for (bl,chan) in self.blChanPairs.keys()]) / 5000.0
        return noiseVariance

    #############################################
    #   Generalized Functions
    #############################################

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.lincalIterations = 0
        self.evWarned = False
        self.blChanPairs = {}
        self.uSamples = []
        self.previousStep = None

    def countBins(self):
        """Counts baselines, channels, visibilities, and uBins for each time these things change."""
        self.bls = list(set([bl for (bl,chan) in self.blChanPairs.keys()]))
        self.chans = sorted(np.unique(np.asarray([chan for (bl,chan) in self.blChanPairs.keys()])))
        self.nChans = len(self.chans)
        self.chanIndices = {chan: i for i,chan in enumerate(self.chans)}
        self.nVis = len(self.blChanPairs)
        try: 
            self.meanAbsVis = np.mean([np.abs(entry['vis']) for entry in self.blChanPairs.values()])
            self.nCoarseBins = len(self.coarseus)
            self.nFineBins = len(self.fineus)
        except: pass

    def setupBaselineChannelPairs(self, freqs, bls, chan2FreqDict, bl2SepDict, redundancyDict):
        """Setup function. Takes a list frequencies (in GHz) and a dictionary that converts channel number to frequecy. 
        Also takes a list of baselines, anda  way to convert them into numpy arrays in meters. Currently u is 1D and not 3D."""
        self.freqs = freqs
        self.chans = sorted(np.unique(np.asarray([chan for chan in chan2FreqDict.keys()])))
        self.chan2FreqDict, self.bl2SepDict = chan2FreqDict, bl2SepDict
        freqChanPairs = sorted([(chan2FreqDict[chan], chan) for chan in self.chans])
        for f,chan in freqChanPairs:
            for bl in bls: self.blChanPairs[(bl,chan)] = {'u': freqs[chan]*bl2SepDict[bl][0], 'bl': bl, 'chan': chan, 'freq': chan2FreqDict[chan]} #TODO: generalize to 3D
        self.countBins()
        self.redundancyDict = redundancyDict
        if self.verbose: print "Examining " + str(self.nVis) + " visibilities from " + str(self.nChans) + " channels..."

    def loadData(self, data, t):
        """Loads data from a data object created from one or more omnical .npzs. Takes t, the integration index into the data."""
        self.t = t
        for (bl,chan),entry in self.blChanPairs.items(): entry['vis'] = data[bl][t][chan]

    def flagChans(self, flaggedChans=[]):
        """Flags channels (and removed corresponding visibilities from blChanPairs) if the data are all 0 or if the channel is passed as flagged."""
        for chan in self.chans:
            if np.all([self.blChanPairs[(bl,chan)]['vis'] == 0.0 for bl in self.bls]): 
                for bl in self.bls: del self.blChanPairs[(bl,chan)]
        for bl,chan in self.blChanPairs.keys(): 
            if chan in flaggedChans: del self.blChanPairs[(bl,chan)]
        self.countBins()
        if self.verbose: print "    After flagging, " + str(self.nVis) + " visibilities from " + str(self.nChans) + " channels remain."

    def setupSkyFreqBasis(self, skyFreqOrder, skyFreqBasisFile):
        """This sets up the basis functions for the frequency-dependence of the underlying FTed sky."""
        self.skyFreqOrder = skyFreqOrder
        if skyFreqBasisFile is not None: self.skyFreqBasis = np.load(skyFreqBasisFile)
        else: self.skyFreqBasis = np.array([np.log(self.freqs/.145)**order for order in range(skyFreqOrder+1)]).T

    def defineuSampling(self, desiredDeltau, uMin, uMax, beamFullRange):
        """This function cuts visibilities with too small or large u, then figures out a list of bed-of-nails u modes that need to be estimated."""
        
        for (bl,chan),entry in self.blChanPairs.items():
            if entry['u'] < uMin or entry['u'] > uMax: del self.blChanPairs[(bl,chan)]

        #TODO: generalize to 2D arrays
        self.fineDu = np.min(np.diff(sorted([self.bl2SepDict[bl][0] for bl in self.bls]))) * (self.freqs[1]-self.freqs[0])
        self.coarseDu = np.floor(desiredDeltau / self.fineDu) * self.fineDu
        self.padding = (np.floor(desiredDeltau / self.fineDu)-1)/2.0

        self.beamFullRange = beamFullRange
        self.beamCoarseBins = int(np.ceil(self.beamFullRange / self.coarseDu))
        self.beamCoarseBins = self.beamCoarseBins + (self.beamCoarseBins+1)%2
        self.beamFineBins = int((self.beamCoarseBins+1) * (2*self.padding + 1))

        for entry in self.blChanPairs.values(): entry['u'] = int(np.round(entry['u'] / self.fineDu)) * self.fineDu
        allus = [entry['u'] for entry in self.blChanPairs.values()]

        coarseBuffer = (self.beamCoarseBins-1)/2 * self.coarseDu
        self.coarseus = np.arange(np.min(allus) - coarseBuffer, np.max(allus) + coarseBuffer + 1e-10, self.coarseDu)
        self.fineus = np.arange(self.coarseus[0] - np.floor(self.padding)*self.fineDu, self.coarseus[-1] + np.ceil(self.padding)*self.fineDu + 1e-10, self.fineDu) 
        self.coarseDeltaus = np.arange(-(self.beamCoarseBins-1)/2*self.coarseDu, (self.beamCoarseBins)/2*self.coarseDu+1e-10, self.coarseDu)
        self.fineDeltaus = np.arange(self.coarseDeltaus[0]-int(np.floor(self.padding))*self.fineDu, self.coarseDeltaus[-1]+int(np.ceil(self.padding))*self.fineDu+1e-10, self.fineDu)

        #Make a dictionary that maps coarse u values to find indices. Also, each entry knows about the index into fineus and the nearest coarseus
        self.coarseu2fineIndex = {coarseBin: int(i * (self.padding*2 + 1) + np.floor(self.padding)) for i, coarseBin in enumerate(self.coarseus)}
        for entry in self.blChanPairs.values(): entry['fineuBin'] = int(np.round((entry['u'] - self.fineus[0])/self.fineDu))
        for entry in self.blChanPairs.values(): entry['coarseuBin'] = int(np.round((entry['u'] - self.coarseus[0])/self.coarseDu))
        self.countBins()
        if self.verbose: 
            print "    After cuts based on u, " + str(self.nVis) + " visibilities from " + str(self.nChans) + " channels remain."
            print "    With " + str(self.coarseDu) + " wavelength spacing, there are " + str(self.nCoarseBins) + " coarse u bins."

    def computeModel(self, blChanPair, factorFunctions, factorParams):
        """This function evaluates the model, the product of the factorFunctions using the factorParams, evaluated for this bl-chan pair."""
        model = 1.0 + 0.0j
        for func,params in zip(factorFunctions, factorParams): model *= func(blChanPair, params)
        return model

    def computeModelsAndErrors(self, factorFunctions, factorParams):
        """This function returns a list of all errors by evaluating the model for all blChanPairs in the .items() order.
        The function also saves the model in each blChanPair as a ['model'] dictionary entry. """
        errors = np.zeros(self.nVis, dtype=complex)
        for n,entry in enumerate(self.blChanPairs.values()):
            entry['model'] = self.computeModel(entry, factorFunctions, factorParams)
            errors[n] = entry['vis'] - entry['model']
        return errors

    def performLogcal(self):
        """Performs logcal and returns normalized first guesses at beta and Sigma (but only to 0th order)."""
        chan2Col = {chan: col for col,chan in zip(range(0,self.nChans),self.chans)}
        uBin2Col = {uBin: col for col,uBin in zip(range(self.nChans,self.nChans+self.nCoarseBins),range(self.nCoarseBins))}
        Acoeffs, rowIndices, colIndices = [np.zeros(self.nVis*2) for i in range(3)] 
        for n,((bl,chan), entry) in enumerate(self.blChanPairs.items()):
                rowIndices[2*n:2*n+2] = n
                colIndices[2*n:2*n+2] = [chan2Col[chan], uBin2Col[entry['coarseuBin']]]
                Acoeffs[2*n:2*n+2] = [1.0, 1.0]
        
        self.logcalA = csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(self.nVis, self.nChans + self.nCoarseBins))
        self.logcalAtA = self.logcalA.conjugate().transpose().dot(self.logcalA).toarray()
        logcalZeroEVs = len(self.logcalAtA) - np.linalg.matrix_rank(self.logcalAtA)
        if logcalZeroEVs > 1: print "    WARNING: Logcal's AtA has " + str(logcalZeroEVs) + " zero eigenvalues. It should have just 1."

        y = np.asarray([np.log(entry['vis']) for entry in self.blChanPairs.values()])
        xhatReal, xhatImag = [np.linalg.pinv(self.logcalAtA).dot((self.logcalA.conjugate().T).dot(thisy)) for thisy in (np.real(y), np.imag(y))]
        return self.initializeLogcalResults(np.exp(xhatReal + 1.0j*xhatImag))

    def generalizedLincalIteration(self, factorFunctions, factorParams, noiseVariance, alpha=1.0):
        """TODO: document this"""
        self.lincalIterations += 1
        
        #Convert parameters to column numbers:
        nCols, toCols = 0, []
        for params in factorParams: 
            toCols.append({index:col for index,col in enumerate(range(nCols, nCols+len(params)))})
            nCols += len(params)

        #Calculate A matrix entries
        coeffs, rowIndices, colIndices = [], [], []
        for row, entry in enumerate(self.blChanPairs.values()):
            model = self.computeModel(entry, factorFunctions, factorParams)
            for func,params,toCol in zip(factorFunctions, factorParams,toCols):
                derivatives, indices = func(entry, params, derivative=True) 
                colIndices.append([toCol[index] for index in indices])
                rowIndices.append([row] * len(indices))
                coeffs.append(np.asarray(derivatives) * model / (func(entry, params)+1e-15))
        rowIndices = np.fromiter(chain.from_iterable(rowIndices), dtype='int')
        colIndices = np.fromiter(chain.from_iterable(colIndices), dtype='int')
        coeffs = np.fromiter(chain.from_iterable(coeffs), dtype='complex128')

        #Generate A, AtNinvA
        self.A = csr_matrix((coeffs,(rowIndices,colIndices)), shape=(self.nVis, nCols))
        Ninv = csr_matrix(((noiseVariance)**-1, (np.arange(self.nVis), np.arange(self.nVis))))
        self.AtNinvA = (self.A.conjugate().transpose().dot(Ninv)).dot(self.A).toarray()
        #TODO: put back in missing eigenvalue warnings
        if True:#self.evWarned or self.lincalIterations == 1:
            zeroEVs = len(self.AtNinvA) - np.linalg.matrix_rank(self.AtNinvA + self.dampingFactor*np.diag(np.diag(self.AtNinvA)))
            if zeroEVs > 1: 
                self.evWarned = True
                print "        WARNING: Lincals's AtNinvA has " + str(zeroEVs) + " zero eigenvalues. It should have just 1."
                print "        CONDITION NUMBER: " + str(np.linalg.cond(self.AtNinvA + self.dampingFactor*np.diag(np.diag(self.AtNinvA))))
        



        errors = self.computeModelsAndErrors(factorFunctions, factorParams)
        currentChiSqPerDoF = np.average(np.abs(errors)**2 / noiseVariance)
        while True:

            M = self.AtNinvA + self.dampingFactor*np.diag(np.diag(self.AtNinvA))
            xHat = np.linalg.pinv(M).dot(self.A.T.conjugate().dot(Ninv.dot(errors)))
            oldParams = copy.deepcopy(factorParams)
            proposedStep = alpha * xHat
            if self.previousStep is not None: 
                stepCosine = np.real(np.dot(proposedStep, self.previousStep.conjugate()))
                stepCosine /= np.linalg.norm(proposedStep) * np.linalg.norm(self.previousStep)
                print 'stepCosine: ', stepCosine
            else: stepCosine = 0.0

            #update params and compute Chi^2
            for params, toCol in zip(factorParams, toCols): 
                for i in range(len(params)): params[i] = params[i] + alpha*xHat[toCol[i]]
            if self.dampingFactor > 1e16: factorParams = oldParams #calls it quits 
            self.renormalize(factorFunctions, factorParams, absMeanGoal = self.meanAbsVis)
            chiSqPerDoF = np.average(np.abs(self.computeModelsAndErrors(factorFunctions, factorParams))**2 / noiseVariance)
            if self.dampingFactor > 1e16: break

            if chiSqPerDoF < currentChiSqPerDoF or (1-stepCosine)**2 * chiSqPerDoF <= currentChiSqPerDoF: #bold acceptance criterion
                if self.dampingFactor > 1e-20: self.dampingFactor /= 5.0#2.0
                self.previousStep = proposedStep
                break
            else: 
                self.dampingFactor *= 2.0#1.5 #as suggested by http://arxiv.org/pdf/1201.5885v1.pdf
                factorParams = oldParams
                print self.dampingFactor, chiSqPerDoF, (1-stepCosine)**2 * chiSqPerDoF - currentChiSqPerDoF, stepCosine
        
        return factorParams, chiSqPerDoF


    def renormalizeNoise(self, factorFunctions, factorParams, noiseVariance):
        """Returns a new noise variance that has been renormalized such that the median observed error reflects the median noise covariance. """
        errors = self.computeModelsAndErrors(factorFunctions, factorParams)
        return noiseVariance * np.median(np.abs(errors)**2) / np.median(noiseVariance)

    def SortedEigensystem(self, matrix):
        """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
        evals,evecs = np.linalg.eig(matrix)
        indices = np.argsort(np.abs(evals))[::-1]   
        return evals[indices], evecs[:,indices]    


    # def performLincalIteration(self, betas, Sigmas, noiseVariance, alpha=1.0):
    #     """Performs lincal, using betas and Sigmas fromt the previous calculation. Now Sigmas is an array of size nuBins x (2*SigmaOrder+1).
    #     Also takes in noise variance model and alpha, which slows down convergence to potentially avoid false minima. 
    #     Returns updated and renormalized betas, Sigmas, and a chi^2 / DoF."""
    #     self.lincalIterations += 1
    #     chan2Col = {chan: col for col,chan in enumerate(self.chans)}
    #     uBin2Cols = {uBin: range(self.nChans+binNum*(2*self.SigmaOrder+1), self.nChans+(binNum+1)*(2*self.SigmaOrder+1)) for binNum,uBin in enumerate(self.uBins)}
    #     nTerms = 1+2*self.SigmaOrder+1
    #     coeffs, rowIndices, colIndices = np.zeros(self.nVis*nTerms,dtype=complex), np.zeros(self.nVis*nTerms), np.zeros(self.nVis*nTerms)
    #     for n,((bl,chan), entry) in enumerate(self.blChanPairs.items()):
    #         uBin = entry['uBin']
    #         rowIndices[nTerms*n:nTerms*(n+1)] = [n]*nTerms
    #         colIndices[nTerms*n:nTerms*(n+1)] = [chan2Col[chan]] + uBin2Cols[uBin]
    #         deltau = entry['u'] - self.uBinCenters[uBin]
    #         firstTerm =  [betas[chan] * self.complexFourier(deltau, Sigmas[uBin])]
    #         coeffs[nTerms*n:nTerms*(n+1)] = np.asarray(firstTerm + [betas[chan] * np.exp(1.0j*np.pi * m * deltau / self.uBinSize) for m in range(-self.SigmaOrder, self.SigmaOrder+1)])

    #     self.A = csr_matrix((coeffs,(rowIndices,colIndices)))
    #     Ninv = csr_matrix(((noiseVariance)**-1, (np.arange(self.nVis), np.arange(self.nVis))))
    #     self.AtNinvA = (self.A.conjugate().transpose().dot(Ninv)).dot(self.A).toarray()
    #     if self.evWarned or self.lincalIterations == 1:
    #         self.evWarned = True
    #         zeroEVs = len(self.AtNinvA) - np.linalg.matrix_rank(self.AtNinvA)
    #         if zeroEVs > 1: print "        WARNING: Lincals's AtNinvA has " + str(zeroEVs) + " zero eigenvalues. It should have just 1."
        
    #     xHat = np.linalg.pinv(self.AtNinvA).dot(self.A.T.conjugate().dot(Ninv.dot(self.computeModelErrors(betas, Sigmas))))
    #     newBetas = {chan: betas[chan]*(1+alpha*xHat[chan2Col[chan]]) for chan in self.chans}
    #     newSigmas = {uBin: Sigmas[uBin] + alpha*np.asarray(xHat[uBin2Cols[uBin]]) for uBin in self.uBins}
    #     betas, Sigmas = self.renormalize(newBetas, newSigmas, absMeanGoal=self.meanAbsVis)
    #     chiSqPerDoF = np.average(np.abs(self.computeModelErrors(betas, Sigmas))**2 / noiseVariance) 
    #     return betas, Sigmas, chiSqPerDoF





#############################################
#   I/O Functionality
#############################################

# def to_npz(resultsFile, dataFiles, allBandpasses, meanBandpass, stdBandpass, unflaggedChans, channelRMSs, overallChannelRMS, bandpassFit=None):
#     """This function saves all useful results for plotting and later combining results from multiple files. Including:
#         -- dataFiles: list of data files that went into this
#         -- allBandpasses: numpy array of bandpasses for each time and frequency. Zeros where there are flags.
#         -- meanBandpass: averaged bandpass over all times (flagged channels not in average)
#         -- stdBandpass: averaged bandpass over all times (flagged channels not in std)
#         -- unflaggedChans: list of list of channels that are not flagged for each time (from uCal.chans)
#         -- channelRMSs: numpy array of data - model RMS for each channel and each time
#         -- overallChannelRMS: RMSs combined over all times
#         -- bandPassfit: lower order fit to the bandpass, evaluated at the channels. Optional. """
#     np.savez(resultsFile, dataFiles=dataFiles, allBandpasses=allBandpasses, meanBandpass=meanBandpass, stdBandpass=stdBandpass, 
#         unflaggedChans=unflaggedChans, channelRMSs=channelRMSs, overallChannelRMS=overallChannelRMS, bandpassFit=bandpassFit)
#     print 'Now saving results to ' + resultsFile + '\n'

# def from_npz(resultsFile):
#     """Loads summary results. See definition of to_npz"""
#     r = np.load(resultsFile)
#     return r['dataFiles'], r['allBandpasses'], r['meanBandpass'], r['stdBandpass'], r['unflaggedChans'], r['channelRMSs'], r['overallChannelRMS'], r['bandpassFit']
    

