import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from copy import deepcopy


class InterferometricArray():
    """Class that takes a list of positions and can calcuate baselines and redundancies."""
    
    def __init__(self, positions=[]):
        self.positions = np.array(positions)
        self.nant = len(positions)
        self.antNames = range(self.nant)
    
    def CalculateUBLs(self, precisionFactor=1000000):
        """Finds the baselines, unique baselines, and related dictionaries for indexing."""
        self.blVectors, self.blNamePairs, self.blIndexPairs = [], [], []
        self.index2name = {i:name for i,name in enumerate(self.antNames)}
        self.name2index = {name:i for i,name in enumerate(self.antNames)}
        self.name2pos = {name: self.positions[self.name2index[name]] for name in self.antNames}
        for index1,ant1 in enumerate(self.antNames):
            for index2,ant2 in zip(range(index1+1,self.nant), self.antNames[index1+1:]):
                delta = np.array([int(np.round(precisionFactor*(self.positions[index1][i] - self.positions[index2][i]))) for i in range(3)])
                if delta[1] > 0 or (delta[1] == 0 and delta[0] > 0): 
                    self.blVectors.append(tuple(delta))
                    self.blNamePairs.append((ant1, ant2))
                    self.blIndexPairs.append((index1, index2))
                else: 
                    self.blVectors.append(tuple(-delta))
                    self.blNamePairs.append((ant2, ant1))
                    self.blIndexPairs.append((index2, index1))
        self.ublDict = {}
        for b in range(len(self.blVectors)):
            if self.ublDict.has_key(self.blVectors[b]): self.ublDict[self.blVectors[b]].append(self.blNamePairs[b])
            else: self.ublDict[self.blVectors[b]] = [self.blNamePairs[b]]
        self.blIndexDict = {antPair: i for i,antPair in enumerate(self.blNamePairs)}
        self.names2ublIndex = {antPair: i for i,antPairs in enumerate(self.ublDict.values()) for antPair in antPairs}
        self.indices2ublIndex = {(self.name2index[antPair[0]],self.name2index[antPair[1]]): 
                                 i for i,antPairs in enumerate(self.ublDict.values()) for antPair in antPairs}
        self.ublVectors = np.array([self.name2pos[antList[0][0]]-self.name2pos[antList[0][1]] for antList in self.ublDict.values()])
        self.ublGroups = [antList for antList in self.ublDict.values()]
        print "With", len(self.positions), "antennas there are", len(self.ublDict.items()), "unique baselines."
        self.nbl, self.nubl = len(self.blNamePairs), len(self.ublVectors)

class HexagonalArray(InterferometricArray):
    """Generates a hexagonal array."""
    def __init__(self, separation, hexNum, verbose=False):
        """Creates a hexagonal array with hexNum antennas per side separated by separation."""
        self.hexNum, self.separation, self.verbose = hexNum, separation, verbose
        positions, self.rowIndices, i = [], [], 0        
        for row in range(hexNum-1,-(hexNum),-1):
            indices = []
            for col in range(2*hexNum-abs(row)-1):
                xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*separation;
                yPos = row*separation*3**.5/2;
                positions.append([xPos, yPos, 0])
                indices.append(i); i+=1
            self.rowIndices.append(indices)
        self.positions = np.array(positions) 
        self.nant = len(self.positions)
        self.antNames = range(self.nant)

class InterferometricArrayFromAA(InterferometricArray):
    """Generates an InterferometricArray object from an aipy antenna array obejct. 
    Supports antenna names that are not the same as antenna indices in the list of positions."""
    def __init__(self, aa, xants=[]):
        self.positions = np.array([pos for i,pos in enumerate(aa.antpos_ideal) if pos[2] != -1.0 and i not in xants])
        self.antNames = [i for i,pos in enumerate(aa.antpos_ideal) if pos[2] != -1.0 and i not in xants]
        self.nant = len(self.positions)


class OmniCalibrator():
    """This class contains method functions for the key steps in Omnical and stores relevant information about the array."""
    
    def __init__(self, array, verbose=True):
        self.a = array
        self.verbose = verbose
        antloc = np.array(self.a.positions); antloc[:,2] = 0
        self.Rgains = np.array([np.append(ai,1) for ai in antloc]); #Get the R matrix. R = [r_i 1], where ri are the positions
        self.Mgains = np.linalg.pinv(self.Rgains.T.dot(self.Rgains)).dot(self.Rgains.T)
        self.Rvis = np.hstack((-self.a.ublVectors, np.zeros((self.a.nubl,1))))
        return
    
    def ComputeErrors(self, obsVis, gainSols, visSols):
        """Computes the difference between the calibration model and the observation."""
        modelObs = np.array([gainSols[ant1] * np.conj(gainSols[ant2]) * visSols[self.a.indices2ublIndex[(ant1,ant2)]] 
                             for (ant1,ant2),obs in zip(self.a.blIndexPairs,obsVis)])
        return obsVis - modelObs

    def PerChannelDegeneracyCorrection(self, gainSols, visSols, degenGains):
        """This function fixes the gain and phase degeneracies on a per-channel basis. Vulnerable to phase wrapping and related problems."""
        newGainSols, newVisSols = gainSols.copy(), visSols.copy()

        #Fix amplitudes
        newGainSols = gainSols * np.exp(1.0j * (np.mean(np.angle(degenGains)) - np.mean(np.angle(newGainSols))))
        newGainSols = gainSols / np.mean(np.abs(gainSols)) * np.mean(np.abs(degenGains))
        newVisSols = visSols * np.mean(np.abs(gainSols))**2 / np.mean(np.abs(degenGains))**2
            
        #Fix phases
        degenRemoved = self.Mgains.dot(np.angle(newGainSols)) - self.Mgains.dot(np.angle(degenGains))
        newGainSols = newGainSols * np.exp(-1.0j * self.Rgains.dot(degenRemoved))
        newVisSols = newVisSols * np.exp(-1.0j * self.Rvis.dot(degenRemoved))
        return newGainSols, newVisSols    
    
    def PerformLogcal(self, obsVis, degenGains):
        """Performs logcal using obsVis and self.a. Degeneracies are fixed by using degenGains, degenVis."""
        Acoeffs, Bcoeffs, rowIndices, colIndices = [np.zeros(self.a.nbl*3) for i in range(4)]
        for n,(ant1,ant2) in enumerate(self.a.blIndexPairs):
            rowIndices[3*n:3*n+3] = n
            colIndices[3*n:3*n+3] = [ant1, ant2, self.a.nant + self.a.indices2ublIndex[(ant1,ant2)]]
            Acoeffs[3*n:3*n+3] = [1.0, 1.0, 1.0]
            Bcoeffs[3*n:3*n+3] = [1.0, -1.0, 1.0]

        logcalA = csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(self.a.nbl, self.a.nant + self.a.nubl))
        logcalB = csr_matrix((Bcoeffs,(rowIndices,colIndices)), shape=(self.a.nbl, self.a.nant + self.a.nubl))
        AtA = (logcalA.conj().T.dot(logcalA)).toarray()
        BtB = (logcalB.conj().T.dot(logcalB)).toarray()
        nZeroEigenvalues = [np.sum(np.abs(np.linalg.eigvals(XtX)<1e-10)) for XtX in [AtA, BtB]]
        if not nZeroEigenvalues == [1,3]: print 'WARNING: Array is not omnical-able.' 

        xReal = (np.linalg.pinv(AtA)).dot(logcalA.conj().T.dot(np.real(np.log(obsVis))))
        xImag = (np.linalg.pinv(BtB)).dot(logcalB.conj().T.dot(np.imag(np.log(obsVis))))
        xHat = np.exp(xReal + 1.0j*xImag)
        logcalGainSols, logcalVisSols = xHat[0:self.a.nant], xHat[self.a.nant:]        
        return self.PerChannelDegeneracyCorrection(logcalGainSols, logcalVisSols, degenGains)
    
    def LincalAMatrix(self, gainSols, visSols, realImagMode=False):
        """Calculates A used for lincal as a compressed sparse row matrix."""
        Acoeffs, rowIndices, colIndices = [np.zeros(self.a.nbl*12) for i in range(3)]
        for n,(ant1,ant2) in enumerate(self.a.blIndexPairs):
            rowIndices[12*n:12*n+6] = 2*n
            rowIndices[12*n+6:12*n+12] = 2*n+1
            ublIndex = self.a.indices2ublIndex[(ant1,ant2)]
            colIndices[12*n:12*n+6] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*self.a.nant+2*ublIndex, 2*self.a.nant+2*ublIndex+1]
            colIndices[12*n+6:12*n+12] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*self.a.nant+2*ublIndex, 2*self.a.nant+2*ublIndex+1]
            if realImagMode: #Compute coefficients for Real/Imag version of lincal
                gi0V0 = gainSols[ant1] * visSols[ublIndex]
                gj0starV0 = np.conj(gainSols[ant2]) * visSols[ublIndex]
                gi0gj0star = gainSols[ant1] * np.conj(gainSols[ant2])
                Acoeffs[12*n:12*n+6] = [np.real(gj0starV0), -np.imag(gj0starV0), np.real(gi0V0), 
                                        np.imag(gi0V0), np.real(gi0gj0star), -np.imag(gi0gj0star)]
                Acoeffs[12*n+6:12*n+12] = [np.imag(gj0starV0), np.real(gj0starV0), np.imag(gi0V0), 
                                           -np.real(gi0V0), np.imag(gi0gj0star), np.real(gi0gj0star)]
            else: #Compute coefficients for Amp/Phase version of lincal
                gi0gj0star = gainSols[ant1] * np.conj(gainSols[ant2])
                gi0gj0starVij0 = gi0gj0star * visSols[ublIndex]
                Acoeffs[12*n:12*n+6] = [np.real(gi0gj0starVij0), -np.imag(gi0gj0starVij0), np.real(gi0gj0starVij0), 
                                        np.imag(gi0gj0starVij0), np.real(gi0gj0star), -np.imag(gi0gj0star)]
                Acoeffs[12*n+6:12*n+12] = [np.imag(gi0gj0starVij0), np.real(gi0gj0starVij0), np.imag(gi0gj0starVij0), 
                                           -np.real(gi0gj0starVij0), np.imag(gi0gj0star), np.real(gi0gj0star)]
        return csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(2*self.a.nbl, 2*self.a.nant + 2*self.a.nubl))
    
    def PerformLincal(self, obsVis, gainStart, visStart, degenGains, realImagMode=False, maxIter=100, convCrit=1e-14, divCrit = 1e14):
        """Performs lincal, either in the amp/phase mode or the real/imag mode. Projects out degeneracies and replaces them."""
        gainSols, visSols = gainStart.copy(), visStart.copy()
        if self.verbose: print '\nNow performing Lincal using the', ('Amp/Phs','Re/Imag')[realImagMode], 'method...'
        startingChiSq = np.mean(np.abs(self.ComputeErrors(obsVis, gainSols, visSols))**2) #TODO: update with noise
        prevChiSq = startingChiSq
    
        for iteration in range(maxIter):            
            #Do all the linear algebra
            A = self.LincalAMatrix(gainSols, visSols, realImagMode=realImagMode)
            AtA = (A.conj().T.dot(A)).toarray()
            error = self.ComputeErrors(obsVis, gainSols, visSols)
            y = np.dstack((np.real(error),np.imag(error))).flatten()
            xHat = np.linalg.pinv(AtA).dot(A.conj().T.dot(y))
            
            #Update solutions
            updates = xHat[0::2] + 1.0j*xHat[1::2]
            if realImagMode: newGainSols = gainSols + updates[0:self.a.nant]
            else: newGainSols = gainSols * (1.0+updates[0:self.a.nant])
            newVisSols = visSols + updates[self.a.nant:]

            #Fix degeneracies if things have gone haywire
            chiSq = np.mean(np.abs(self.ComputeErrors(obsVis, newGainSols, newVisSols))**2)
            if chiSq > startingChiSq:
                newGainSols, newVisSols = self.PerChannelDegeneracyCorrection(newGainSols, newVisSols, degenGains)
            convergence = np.linalg.norm(np.append(newGainSols-gainSols,newVisSols-visSols)) / np.linalg.norm(np.append(newGainSols,newVisSols))
            #convergence = np.abs(prevChiSq - chiSq) / chiSq
            if self.verbose: print iteration, ' -- chiSq =', chiSq, '    convCrit =', convergence
            
            #Check for convergence
            if chiSq/startingChiSq > divCrit or iteration == maxIter-1:
                print 'WARNING: Lincal in ' + ('Amp/Phs','Re/Imag')[realImagMode] + ' Mode did not converge in', iteration+1, 'iterations!'
                print '    Stopped at chiSq =', chiSq, '    convCrit =', convergence
                break
            if convergence < convCrit : 
                if self.verbose: print ('Amp/Phs','Re/Imag')[realImagMode] + ' Mode converged in', iteration+1, 'iterations.'
                if self.verbose: print '    Stopped at chiSq =', chiSq, '    convCrit =', convergence
                break
            gainSols, visSols, prevChiSq = newGainSols, newVisSols, chiSq
        
        return self.PerChannelDegeneracyCorrection(newGainSols, newVisSols, degenGains)
    
    def OverallBandpassDegeneracyProjection(self, allGainSols, allVisSols, allGuessGains):
        """This function corrects degeneracies, but only after unwrapping across all channels. """
        finalGainSols, finalVisSols, allDegenGains = [np.array(this) for this in [allGainSols, allVisSols, allGuessGains]]
        #overall phase unwrapping
        phaseShift = np.mean(np.unwrap(np.angle(finalGainSols).T).T,axis=1)
        for ant in range(self.a.nant): finalGainSols[:,ant] *= np.exp(-1.0j * phaseShift)

        allDegens = []
        for chan in range(len(allGainSols)):
            unwrappedGainPhases = np.array([np.unwrap(np.angle(finalGainSols[0:chan+1,ant]))[chan] for ant in range(self.a.nant)])
            unwrappedDegenPhases = np.array([np.unwrap(np.angle(np.array(allDegenGains)[0:chan+1,ant]))[chan] for ant in range(self.a.nant)])
            degenRemoved = self.Mgains.dot(unwrappedGainPhases) - self.Mgains.dot(unwrappedDegenPhases)
            finalGainSols[chan,:] = finalGainSols[chan,:] * np.exp(-1.0j * self.Rgains.dot(degenRemoved))
            finalVisSols[chan,:] = finalVisSols[chan,:] * np.exp(-1.0j * self.Rvis.dot(degenRemoved))

        return finalGainSols, finalVisSols

    def determineConjugationAndExemplars(self, data, v, pol='x'):
        """This function examines the data and guess visibility objects to build up relevant objects that say how the data is stored."""
        self.dataKeyConjugated = {}
        for (ant1,ant2) in self.a.blNamePairs:
            if data.has_key((ant1,ant2)): self.dataKeyConjugated[(ant1,ant2)] = False
            elif data.has_key((ant2,ant1)): self.dataKeyConjugated[(ant1,ant2)] = True
            else: raise KeyError('Cannot find baseline' + str((ant1,ant2)) + " in observed visibility data.")
        
        #Now examine the input visibility quantities
        self.visKey2ublIndex, self.visKeyConjugated = {}, {}
        self.ublIndex2VisKey = [None for i in range(self.a.nubl)]
        for n,ublGroup in enumerate(self.a.ublGroups):
            for (ant1,ant2) in ublGroup:
                if v[pol+pol].has_key((ant1,ant2)): 
                    self.visKey2ublIndex[(ant1,ant2)] = n
                    self.visKeyConjugated[(ant1,ant2)] = False
                    self.ublIndex2VisKey[n] = (ant1,ant2)
                elif v[pol+pol].has_key((ant2,ant1)): 
                    self.visKey2ublIndex[(ant2,ant1)] = n
                    self.visKeyConjugated[(ant2,ant1)] = True
                    self.ublIndex2VisKey[n] = (ant2,ant1)

    def pullOutSingleTime(self, time, data, g, v, pol='x', defaultVis=1.0):
        """This function pulls out observations, gains, and visibilities from the omnical data format. 
        If a particular unique baseline doesn't appear in v, defaultVis is used instead."""
        nint, nchan = data.values()[0].values()[0].shape
        allObsVis = np.ones((nchan, self.a.nbl), dtype=complex)
        allGains = np.ones((nchan, self.a.nant), dtype=complex)
        allVis = np.ones((nchan, self.a.nubl), dtype=complex)

        #Extract observed visbility data
        for n,(ant1,ant2) in enumerate(self.a.blNamePairs):
            if self.dataKeyConjugated[(ant1,ant2)]: allObsVis[:,n] =  np.conj(data[(ant2,ant1)][pol+pol][time])
            else: allObsVis[:,n] =  data[(ant1,ant2)][pol+pol][time]
        
        #Extract gains
        for n,ant in enumerate(self.a.antNames): allGains[:,n] = g[pol][ant][time]
        
        #Extract unique visibilities
        for n in range(self.a.nubl):
            visKey = self.ublIndex2VisKey[n]
            if visKey is None: allVis[:,n] *= defaultVis
            elif self.visKeyConjugated[visKey]: allVis[:,n] = np.conj(v[pol+pol][visKey][time])
            else: allVis[:,n] = v[pol+pol][visKey][time]

        return allObsVis, allGains, allVis
    
    def putBackSingleTime(self, time, chans, gSol, vSol, allGainSols, allVisSols, pol='x'):
        """For the single time slice and the specified channels, put gain/vis solutions into omnical format."""
        for n,ant in enumerate(self.a.antNames):
            gSol[pol][ant][time][chans] = allGainSols[:,n]

        for n in range(self.a.nubl):
            visKey = self.ublIndex2VisKey[n]
            if visKey is not None:
                if self.visKeyConjugated[visKey]: vSol[pol+pol][visKey][time][chans] = np.conj(allVisSols[:,n])
                else: vSol[pol+pol][visKey][time][chans] = allVisSols[:,n]
    
    def Omnical(self, data, g0, v0, pol='x', times=[], chans=[], flags=None):
        """Takes in data, g0, v0 in the Omnical dictionary format. Degeneracies set by g0, v0. 
        Only works for xx or yy polarization.
        First-Cal solution should be applied to the data before calibrations starts.
        To calibrate only a subset of channels or times, specify their indices as a list. 
        The uncalibrated portion will be taken from g0 and v0."""
        
        gSol, vSol = deepcopy(g0), deepcopy(v0)
        nint, nchan = data.values()[0].values()[0].shape
        if len(times) == 0: times = range(nint)
        if len(chans) == 0: chans = range(nchan)
        self.determineConjugationAndExemplars(data, v0, pol=pol)
        
        for time in times:            
            allGainSols = np.zeros((len(chans),self.a.nant),dtype=complex)
            allVisSols = np.zeros((len(chans),self.a.nubl),dtype=complex)
            allObsVis, allGuessGains, allGuessVis = self.pullOutSingleTime(time, data, g0, v0, pol=pol)
            
            for n,chan in enumerate(chans):
                obsVis, guessGains = allObsVis[chan], allGuessGains[chan]
                if n==0: 
                    gainStart, visStart = cal.PerformLogcal(obsVis, guessGains)        
                else: gainStart, visStart = allGainSols[n-1], allVisSols[n-1]
                gainSols, visSols = cal.PerformLincal(obsVis, gainStart, visStart, guessGains, realImagMode=True, maxIter=100)
                allGainSols[n,:], allVisSols[n,:] = gainSols, visSols 
                #print np.mean(np.abs(cal.ComputeErrors(obsVis, gainSols, visSols))**2)

            allGainSols, allVisSols = cal.OverallBandpassDegeneracyProjection(allGainSols, allVisSols, allGuessGains)
            self.putBackSingleTime(time, chans, gSol, vSol, allGainSols, allVisSols, pol=pol)
            
        return gSol, vSol