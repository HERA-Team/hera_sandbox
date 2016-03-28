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

#############################################
#   Set-up script specific to PAPER 128
#############################################

pol='xx'
seps = np.arange(1,16) #* 15m  = baseline lengths
dx = 15.0
separations = dx*seps
#freqs = np.arange(.1,.2,.1/203)
freqs = np.arange(.1,.2,.1/203)[160:]
chans = range(len(freqs))
chan2FreqDict = {chan: freqs[chan] for chan in chans}
uTolerance = 15.0 * 15 * (freqs[1]-freqs[0]) / 3e8 * 1e9  

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,dataFiles = o.parse_args(sys.argv[1:])

if len(sys.argv[1:]) == 0: #This option is for testing
    #dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")] 
    dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.57058.xx.uvcRRE.npz")] 
 
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


#############################################
#   uCal Core Classes
#############################################

#TODO needs testing on a 2D or 3D array. In particular, there's a question of whether the code that skips entries to make the loop fasterever misses a good pair
class uCalReds():
    """This class takes a list of baselines and channels, and a way to convert them into physical quantities (freqs in GHz, 
    baselines as numpy arrays in meters), and saves a dictionary that maps (ch1,bl1,ch2,bl2) to u and deltau. 
    Physical separations must all have the same orentiation, otherwise this won't work.
    Only includes baselines-frequency pairs obeying the Deltau threshold (default 0.3)."""
    def __init__(self, bls, chans, bl2SepDict, chan2FreqDict, maxDeltau = .3):
        print 'Now finding all baseline/frequency pairs...'
        self.maxDeltau = maxDeltau
        self.blChanPairs = {}
        freqChanPairs = sorted([[chan2FreqDict[chan], chan] for chan in chans])
        sepBLPairs = sorted([[np.asarray(bl2SepDict[bl]), bl] for bl in bls], key=lambda pair: np.linalg.norm(pair[0])) #sorts by length
        for i,(f1,ch1) in enumerate(freqChanPairs):
            for f2,ch2 in freqChanPairs[i+1:]:
                for j,(sep2,bl2) in enumerate(sepBLPairs):
                    for sep1,bl1 in sepBLPairs[j+1:]:
                        deltau = np.linalg.norm(f1*sep1 - f2*sep2) * 1.0e9 / scipy.constants.c
                        if deltau < maxDeltau and not self.blChanPairs.has_key((ch1,bl1,ch2,bl2)) and not self.blChanPairs.has_key((ch2,bl2,ch1,bl1)):
                            u = (f1*sep1 + f2*sep2)/2.0 * 1e9/ scipy.constants.c
                            self.blChanPairs[(ch1,bl1,ch2,bl2)] = (u, deltau)
        print "   " + str(len(self.blChanPairs)) + " baseline/frequency pairs identified with delta u < " + str(maxDeltau)

    def applyuCut(self, uMin=25, uMax=150):
        for key,value in self.blChanPairs.items():
            if np.linalg.norm(value[0]) < uMin or np.linalg.norm(value[0]) > uMax: del[self.blChanPairs[key]]
        print "   " + str(len(self.blChanPairs)) + " baseline/frequency pairs remain after requiring " + str(uMin) + " < u < " + str(uMax)

    def applyChannelFlagCut(self, flaggedChannels):
        print 'applyChannelFlagCut has not been implemented!!!'



class uCalibrator():
    """docstring for uCalibrator()"""
    def __init__(self, uReds, uBinSize = .72**.5, duBinSize = 5.0/203):
        #Internal format for blChanPairs with relevant info about binning
        self.blChanPairs = {key: {'u': u, 'du': du} for key,(u,du) in uReds.blChanPairs.items()}
        
        #Determine binning: first assign integers to u and du
        allus = np.asarray([value['u'] for key,value in self.blChanPairs.items()])
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
        self.duBinCenters = [np.mean(self.duBinCenters[duBin]) for duBin in self.duBins]
        
        #Now create logcal matrices



#############################################
#   uCal Script
#############################################


uReds = uCalReds(bls, chans, bl2SepDict, chan2FreqDict, maxDeltau=.3)
uReds.applyuCut(uMin=25, uMax=150)
#uReds.applyChannelFlagCut(flaggedChannels) 


uCal = uCalibrator(uReds)

