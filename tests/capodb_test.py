#!/usr/bin/env python
#
#  db_test.py
#  
#
#  Created by Danny Jacobs on 3/12/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import unittest
import capodb as cdb, aipy as a, numpy as n

class TestCapodb(unittest.TestCase):
    def setUp(self):
        lat, lon = '45:00', '90:00'
        freqs = a.loc.get_freqs(0.1,0.1/256,256)
        beam = a.fit.BeamFlat(freqs)
        self.aa = cdb.AntennaArray((lat,lon), [],beam=beam)
    def testdbinit(self):
        """Testing initialization from paper_test1.db"""
        prms = self.aa.get_params()
        self.aa[0].set_params(prms[0])
        self.assertEqual(prms[0]['x'], 8.98)
    def testsavereload(self):
        """Testing save/reload: """
        oldamp = self.aa[0].amp
        self.aa[0].amp += 1
        self.aa[0].save(self.aa.c,self.aa.aa_db_params)
        self.aa[0].prm_update(self.aa.c,self.aa.aa_db_params)
        prms2 = self.aa.get_params()
        self.assertAlmostEqual(self.aa[0].amp, oldamp+1, 10)

if __name__ == '__main__':
    unittest.main()
