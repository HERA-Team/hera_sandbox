import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

class InterferometricArray():
    """Class that takes a list of positions and can calcuate baselines and redundancies."""
    
    def __init__(self, positions=[]):
        self.positions = np.array(positions)
        self.nant = len(positions)
    
    def CalculateUBLs(self, precisionFactor=1000000):
        """Finds the baselines, unique baselines, and related dictionaries for indexing."""
        self.blVectors, self.blPairs = [], []
        for ant1 in range(self.nant):
            for ant2 in range(ant1+1,self.nant):
                delta = np.array([int(np.round(precisionFactor*(self.positions[ant1][i] - self.positions[ant2][i]))) for i in range(3)])
                if delta[1] > 0 or (delta[1] == 0 and delta[0] > 0): 
                    self.blVectors.append(tuple(delta))
                    self.blPairs.append((ant1, ant2))
                else: 
                    self.blVectors.append(tuple(-delta))
                    self.blPairs.append((ant2, ant1))
        self.ublDict = {}
        for b in range(len(self.blVectors)):
            if self.ublDict.has_key(self.blVectors[b]): self.ublDict[self.blVectors[b]].append(self.blPairs[b])
            else: self.ublDict[self.blVectors[b]] = [self.blPairs[b]]
        self.blIndexDict = {antPair: i for i,antPair in enumerate(self.blPairs)}
        self.ublIndexDict = {antPair: i for i,antPairs in enumerate(self.ublDict.values()) for antPair in antPairs }
        self.ublVectors = np.array([self.positions[antList[0][0]]-self.positions[antList[0][1]] for antList in self.ublDict.values()])
        self.ublGroups = [antList for antList in self.ublDict.values()]
        print "With", len(self.positions), "antennas there are", len(self.ublDict.items()), "unique baselines."
        self.nbl, self.nubl = len(self.blPairs), len(self.ublVectors)

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

class OmniCalibrator():
    """This class contains method functions for the key steps in Omnical and stores relevant information about the array."""
    
    def __init__(self, array, verbose=True):
        self.a = array
        self.verbose = verbose
        antloc = np.array(self.a.positions); antloc[:,2] = 0
        R = np.array([np.append(ai,1) for ai in antloc]); #Get the R matrix. R = [r_i 1], where ri are the positions
        self.R2 = np.vstack((R, np.hstack((-self.a.ublVectors, np.zeros((self.a.nubl,1))))))
        self.M2 = np.linalg.pinv(self.R2.T.dot(self.R2)).dot(self.R2.T)
        return
    
    def ComputeErrors(self, obsVis, gainSols, visSols):
        """Computes the difference between the calibration model and the observation."""
        modelObs = np.array([gainSols[ant1] * np.conj(gainSols[ant2]) * visSols[self.a.ublIndexDict[(ant1,ant2)]] 
                             for (ant1,ant2),obs in zip(self.a.blPairs,obsVis)])
        return obsVis - modelObs

    def PerChannelDegeneracyCorrection(self, gainSols, visSols, degenGains, degenVis):
        """This function fixes the gain and phase degeneracies on a per-channel basis. Vulnerable to phase wrapping and related problems."""
        newGainSols, newVisSols = gainSols.copy(), visSols.copy()

        #Fix amplitudes
        newGainSols = gainSols * np.exp(1.0j * (np.mean(np.angle(degenGains)) - np.mean(np.angle(newGainSols))))
        newGainSols = gainSols / np.mean(np.abs(gainSols)) * np.mean(np.abs(degenGains))
        newVisSols = visSols * np.mean(np.abs(gainSols))**2 / np.mean(np.abs(degenGains))**2
        
        #Fix phases
        fullPhaseDegen = self.R2.dot(self.M2)
        fullPhaseProj = np.eye(self.a.nant + self.a.nubl) - fullPhaseDegen
        newPhases = fullPhaseProj.dot(np.angle(np.append(newGainSols, newVisSols))) + fullPhaseDegen.dot(np.angle(np.append(degenGains,degenVis)))
        newSols = np.abs(np.append(newGainSols, newVisSols)) * np.exp(1.0j * newPhases)
        newGainSols, newVisSols = newSols[0:self.a.nant], newSols[self.a.nant:]
        return newGainSols, newVisSols
    
    def PerformLogcal(self, obsVis, degenGains, degenVis):
        """Performs logcal using obsVis and self.a. Degeneracies are fixed by using degenGains, degenVis."""
        Acoeffs, Bcoeffs, rowIndices, colIndices = [np.zeros(self.a.nbl*3) for i in range(4)]
        for n,(ant1,ant2) in enumerate(self.a.blPairs):
            rowIndices[3*n:3*n+3] = n
            colIndices[3*n:3*n+3] = [ant1, ant2, self.a.nant + self.a.ublIndexDict[(ant1,ant2)]]
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
        return self.PerChannelDegeneracyCorrection(logcalGainSols, logcalVisSols, degenGains, degenVis)
    
    def LincalAMatrix(self, gainSols, visSols, realImagMode=False):
        """Calculates A used for lincal as a compressed sparse row matrix."""
        Acoeffs, rowIndices, colIndices = [np.zeros(self.a.nbl*12) for i in range(3)]
        for n,(ant1,ant2) in enumerate(self.a.blPairs):
            rowIndices[12*n:12*n+6] = 2*n
            rowIndices[12*n+6:12*n+12] = 2*n+1
            ublIndex = self.a.ublIndexDict[(ant1,ant2)]
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
    

    def PerformLincal(self, obsVis, gainStart, visStart, degenGains, degenVis, realImagMode=False, maxIter=100, convCrit=1e-14, divCrit = 1e14):
        """Performs lincal, either in the amp/phase mode or the real/imag mode. Projects out degeneracies and replaces them."""
        gainSols, visSols = gainStart.copy(), visStart.copy()
        if self.verbose: print '\nNow performing Lincal using the', ('Amp/Phs','Re/Imag')[realImagMode], 'method...'
        startingChiSq = np.mean(np.abs(self.ComputeErrors(obsVis, gainSols, visSols))**2) #TODO: update with noise
    
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
                newGainSols, newVisSols = self.PerChannelDegeneracyCorrection(newGainSols, newVisSols, degenGains, degenVis)
            convergence = np.linalg.norm(np.append(newGainSols-gainSols,newVisSols-visSols)) / np.linalg.norm(np.append(newGainSols,newVisSols))
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
            gainSols, visSols = newGainSols, newVisSols
        
        return self.PerChannelDegeneracyCorrection(newGainSols, newVisSols, degenGains, degenVis)
    
    def OverallBandpassDegeneracyProjection(self, allGainSols, allVisSols, allGuessGains, allGuessVis):
        """This function corrects degeneracies, but only after unwrapping across all channels. """
        finalGainSols, finalVisSols, allDegenGains, allDegenVis = [np.array(this) for this in [allGainSols, allVisSols, allGuessGains, allGuessVis]]
        #overall phase unwrapping
        phaseShift = np.mean(np.unwrap(np.angle(finalGainSols).T).T,axis=1)
        for ant in range(self.a.nant): finalGainSols[:,ant] *= np.exp(-1.0j * phaseShift)

        #tip/tilt/phase unwrapping
        fullPhaseDegen = self.R2.dot(self.M2)
        fullPhaseProj = np.eye(self.a.nant + self.a.nubl) - fullPhaseDegen
        allDegens = []
        for chan in range(len(allGainSols)):
            unwrappedGainSols = np.array([np.unwrap(np.angle(finalGainSols[0:chan+1,ant]))[chan] for ant in range(self.a.nant)])
            unwrappedVisSols = np.array([np.unwrap(np.angle(finalVisSols[0:chan+1,ubl]))[chan] for ubl in range(self.a.nubl)])
            unwrappedDegenGains = np.array([np.unwrap(np.angle(np.array(allDegenGains)[0:chan+1,ant]))[chan] for ant in range(self.a.nant)])
            unwrappedDegenVis = np.array([np.unwrap(np.angle(np.array(allDegenVis)[0:chan+1,ubl]))[chan] for ubl in range(self.a.nubl)])
            
            newPhases = fullPhaseProj.dot(np.append(unwrappedGainSols, unwrappedVisSols))
            newPhases += fullPhaseDegen.dot(np.append(unwrappedDegenGains,unwrappedDegenVis))
            newSols = np.abs(np.append(finalGainSols[chan,:], finalVisSols[chan,:])) * np.exp(1.0j * newPhases)
            newGainSols, newVisSols = newSols[0:self.a.nant], newSols[self.a.nant:]
            finalGainSols[chan,:], finalVisSols[chan,:] = newGainSols, newVisSols
            
            #Figure out the degenerate part for everyone
            unwrappedGainSols = np.array([np.unwrap(np.angle(finalGainSols[0:chan+1,ant]))[chan] for ant in range(self.a.nant)])
            unwrappedVisSols = np.array([np.unwrap(np.angle(finalVisSols[0:chan+1,ubl]))[chan] for ubl in range(self.a.nubl)])
            allDegens.append(self.M2.dot(np.append(unwrappedGainSols, unwrappedVisSols)))

        #fix channel 0
        degeneratePart =  self.R2.dot(np.median(allDegens,axis=0))
        newPhases = fullPhaseProj.dot(np.append(np.angle(finalGainSols[0,:]), np.angle(finalVisSols[0,:])))
        newPhases += self.R2.dot(np.median(allDegens,axis=0))
        newSols = np.abs(np.append(finalGainSols[0,:], finalVisSols[0,:])) * np.exp(1.0j * newPhases)
        newGainSols, newVisSols = newSols[0:self.a.nant], newSols[self.a.nant:]
        finalGainSols[0,:], finalVisSols[0,:] = newGainSols, newVisSols
        return list(finalGainSols), list(finalVisSols)