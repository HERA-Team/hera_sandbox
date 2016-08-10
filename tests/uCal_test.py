import unittest
import aipy as a, numpy as np
import capo.uCal as uc
import os

class TestuCalReds(unittest.TestCase):
    def setUp(self):
        pass
        #self.freqs = np.arange(.1,.2,.1/203)
        #self.chans = range(len(self.freqs))
        #self.chan2FreqDict = {chan: self.freqs[chan] for chan in chans}

    def test_init(self):
        freqs = [.1, .2, .3, .4]
        chans = range(len(freqs))
        chan2FreqDict = {chan: freqs[chan] for chan in chans}
        bl2SepDict = {(0,1): np.asarray([45.0,0.0]), (0,2): np.asarray([90.0,0.0])}
        reds = uc.uCalReds(freqs, bl2SepDict.keys(), chan2FreqDict, bl2SepDict)
        self.assertEqual(reds.blChanPairs[0,(0,2),1,(0,1)][0][0], 9.0)
        self.assertEqual(reds.blChanPairs[1,(0,2),3,(0,1)][0][0], 18.0)
        self.assertEqual(reds.blChanPairs[0,(0,2),1,(0,1)][0][1], 0)
        self.assertEqual(reds.blChanPairs[1,(0,2),3,(0,1)][0][1], 0)
        self.assertEqual(len(reds.blChanPairs), 2)

    def test_applyuCut(self):
        freqs = [.1, .2, .3, .4]
        chans = range(len(freqs))
        chan2FreqDict = {chan: freqs[chan] for chan in chans}
        bl2SepDict = {(0,1): np.asarray([45.0,0.0]), (0,2): np.asarray([90.0,0.0])}
        reds = uc.uCalReds(freqs, bl2SepDict.keys(), chan2FreqDict, bl2SepDict)
        reds.applyuCut(uMin=10.0)
        self.assertEqual(len(reds.blChanPairs), 1)
        reds.applyuCut(uMax=1.0)
        self.assertEqual(len(reds.blChanPairs), 0)

    def test_applyChannelFlagCut(self):
        freqs = [.1, .2, .3, .4]
        chans = range(len(freqs))
        chan2FreqDict = {chan: freqs[chan] for chan in chans}
        bl2SepDict = {(0,1): np.asarray([45.0,0.0]), (0,2): np.asarray([90.0,0.0])}
        reds = uc.uCalReds(freqs, bl2SepDict.keys(), chan2FreqDict, bl2SepDict)
        reds.applyChannelFlagCut([0])
        self.assertEqual(len(reds.blChanPairs), 1)
        reds.applyChannelFlagCut([1])
        self.assertEqual(len(reds.blChanPairs), 0)

class TestuCalibrator(unittest.TestCase):
    def setUp(self):
        pass

    def testExactCase(self):
        pass
        #TODO: write a version where we create a set of bl-chan pairs, populate the correlations with exact solutions, and then solve

class TestImportExport(unittest.TestCase):
    def setUp(self):
        uc.save2npz('uCal_unittesting_savetest', 1, 2, 3, 4, 5, 6, 7, 8)
        

    def testImport(self):
        uc.loadAndCombine(['uCal_unittesting_savetest.npz'])

    def tearDown(self):
        os.remove('uCal_unittesting_savetest.npz')


if __name__ == '__main__':
    unittest.main()
