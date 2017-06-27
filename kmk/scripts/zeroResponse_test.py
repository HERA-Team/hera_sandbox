import unittest, numpy as np, aipy
import zeroResponse as zR

class TestMethod(unittest.TestCase):
    def test_makeFlatMap(self):
        Tsky = 1.
        fq = .150
        m = zR.makeFlatMap(64, fq, Tsky=Tsky)
        self.assertTrue(np.all(m.map == m[0])) # check uniform
        px_area = 4*np.pi / m.npix()
        lam = aipy.const.c / (fq * 1e9)
        jy_px = 2*aipy.const.k * Tsky / lam**2 * px_area / 1e-23
        self.assertEqual(m[0], jy_px) # check scaling
    def test_makeSynchMap(self):
        fq = .150
        m = zR.makeSynchMap(64, freq=fq)
        self.assertTrue(np.all(m.map == m[0])) # check uniform
        fq2 = .2
        scalar = (fq2/fq)**-2.5
        jy_scalar = (fq2/fq)**2
        m2 = zR.makeSynchMap(64, freq=fq2)
        self.assertAlmostEqual(m[0] * scalar * jy_scalar, m2[0]) # check scaling in fq

if __name__ == '__main__':
    unittest.main()
