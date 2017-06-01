import capo.redcal as om
import capo, numpy as np
import unittest

np.random.seed(0)

def build_reds(nants):
    reds = {}
    for i in xrange(nants):
        for j in xrange(i+1,nants):
            reds[j-i] = reds.get(j-i,[]) + [(i,j)]
    return [reds[i] for i in reds if len(reds[i]) >= 2]

class TestMethods(unittest.TestCase):
    def test_sim_red_data(self):
        reds = build_reds(10)
        pols = ['xx']
        gains, data = om.sim_red_data(reds, pols, stokes_v=True)
        self.assertEqual(len(gains), 10)
        self.assertEqual(len(data), 45-1)
        for bls in reds:
            bl0 = bls[0]
            ai,aj = bl0
            ans0 = data[bl0+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
            for bl in bls[1:]:
                ai,aj = bl
                ans = data[bl+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
                np.testing.assert_almost_equal(ans0, ans, 7)
        pols = ['xx','yy','xy','yx']
        gains, data = om.sim_red_data(reds, pols, stokes_v=True)
        self.assertEqual(len(gains), 20)
        self.assertEqual(len(data), 4*(45-1))
        for bls in reds:
            bl0 = bls[0]
            ai,aj = bl0
            ans0xx = data[bl0+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
            ans0xy = data[bl0+('xy',)] / (gains[(ai,'x')] * gains[(aj,'y')].conj())
            ans0yx = data[bl0+('yx',)] / (gains[(ai,'y')] * gains[(aj,'x')].conj())
            ans0yy = data[bl0+('yy',)] / (gains[(ai,'y')] * gains[(aj,'y')].conj())
            for bl in bls[1:]:
                ai,aj = bl
                ans_xx = data[bl+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
                ans_xy = data[bl+('xy',)] / (gains[(ai,'x')] * gains[(aj,'y')].conj())
                ans_yx = data[bl+('yx',)] / (gains[(ai,'y')] * gains[(aj,'x')].conj())
                ans_yy = data[bl+('yy',)] / (gains[(ai,'y')] * gains[(aj,'y')].conj())
                np.testing.assert_almost_equal(ans0xx, ans_xx, 7)
                np.testing.assert_almost_equal(ans0xy, ans_xy, 7)
                np.testing.assert_almost_equal(ans0yx, ans_yx, 7)
                np.testing.assert_almost_equal(ans0yy, ans_yy, 7)
        gains, data = om.sim_red_data(reds, pols, stokes_v=False)
        self.assertEqual(len(gains), 20)
        self.assertEqual(len(data), 4*(45-1))
        for bls in reds:
            bl0 = bls[0]
            ai,aj = bl0
            ans0xx = data[bl0+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
            ans0xy = data[bl0+('xy',)] / (gains[(ai,'x')] * gains[(aj,'y')].conj())
            ans0yx = data[bl0+('yx',)] / (gains[(ai,'y')] * gains[(aj,'x')].conj())
            ans0yy = data[bl0+('yy',)] / (gains[(ai,'y')] * gains[(aj,'y')].conj())
            np.testing.assert_almost_equal(ans0xy, ans0yx, 7)
            for bl in bls[1:]:
                ai,aj = bl
                ans_xx = data[bl+('xx',)] / (gains[(ai,'x')] * gains[(aj,'x')].conj())
                ans_xy = data[bl+('xy',)] / (gains[(ai,'x')] * gains[(aj,'y')].conj())
                ans_yx = data[bl+('yx',)] / (gains[(ai,'y')] * gains[(aj,'x')].conj())
                ans_yy = data[bl+('yy',)] / (gains[(ai,'y')] * gains[(aj,'y')].conj())
                np.testing.assert_almost_equal(ans0xx, ans_xx, 7)
                np.testing.assert_almost_equal(ans0xy, ans_xy, 7)
                np.testing.assert_almost_equal(ans0yx, ans_yx, 7)
                np.testing.assert_almost_equal(ans0yy, ans_yy, 7)

class TestRedundantInfo(unittest.TestCase):
    def test_build_eq(self):
        reds = build_reds(3)
        bls = reduce(lambda x,y: x+y, reds)
        info = om.RedundantInfo(reds)
        eqs = info.build_eqs(bls, ['xx'], stokes_v=True)
        self.assertEqual(len(eqs), 2)
        self.assertEqual(eqs['g0x * g1x_ * u0xx'], (0,1,'xx'))
        self.assertEqual(eqs['g1x * g2x_ * u0xx'], (1,2,'xx'))
        pols = ['xx','yy','xy','yx']
        info = om.RedundantInfo(reds)
        eqs = info.build_eqs(bls, pols, stokes_v=True)
        self.assertEqual(len(eqs), 2*4)
        self.assertEqual(eqs['g0x * g1y_ * u0xy'], (0,1,'xy'))
        self.assertEqual(eqs['g1x * g2y_ * u0xy'], (1,2,'xy'))
        self.assertEqual(eqs['g0y * g1x_ * u0yx'], (0,1,'yx'))
        self.assertEqual(eqs['g1y * g2x_ * u0yx'], (1,2,'yx'))
        eqs = info.build_eqs(bls, pols, stokes_v=False)
        self.assertEqual(len(eqs), 2*4)
        self.assertEqual(eqs['g0x * g1y_ * u0xy'], (0,1,'xy'))
        self.assertEqual(eqs['g1x * g2y_ * u0xy'], (1,2,'xy'))
        self.assertEqual(eqs['g0y * g1x_ * u0xy'], (0,1,'yx'))
        self.assertEqual(eqs['g1y * g2x_ * u0xy'], (1,2,'yx'))
    def test_solver(self):
        reds = build_reds(3)
        info = om.RedundantInfo(reds)
        gains,d = om.sim_red_data(reds, ['xx'])
        w = {}
        w = dict([(k,1.) for k in d.keys()])
        def solver(data, wgts, sparse, **kwargs):
            np.testing.assert_equal(data['g0x * g1x_ * u0xx'], d[0,1,'xx'])
            np.testing.assert_equal(data['g1x * g2x_ * u0xx'], d[1,2,'xx'])
            if len(wgts) == 0: return
            np.testing.assert_equal(wgts['g0x * g1x_ * u0xx'], w[0,1,'xx'])
            np.testing.assert_equal(wgts['g1x * g2x_ * u0xx'], w[1,2,'xx'])
            return
        info._solver(solver, d)
        info._solver(solver, d, w)
    def test_logcal(self):
        NANTS = 18
        reds = build_reds(NANTS)
        info = om.RedundantInfo(reds)
        gains,d = om.sim_red_data(reds, ['xx'], gain_scatter=.55)
        w = dict([(k,1.) for k in d.keys()])
        sol = info.logcal(d)
        for i in xrange(NANTS):
            self.assertEqual(sol[(i,'x')].shape, (10,10))
        for bls in reds:
            ubl = sol[bls[0]+('xx',)]
            self.assertEqual(ubl.shape, (10,10))
            for bl in bls:
                d_bl = d[bl+('xx',)]
                mdl = sol[(bl[0],'x')] * sol[(bl[1],'x')].conj() * ubl
                np.testing.assert_almost_equal(np.abs(d_bl), np.abs(mdl), 10)
                np.testing.assert_almost_equal(np.angle(d_bl*mdl.conj()), 0, 10)
    def test_lincal(self):
        NANTS = 18
        reds = build_reds(NANTS)
        info = om.RedundantInfo(reds)
        #gains,d = om.sim_red_data(reds, ['xx'], gain_scatter=.01) # XXX causes svd error
        gains,d = om.sim_red_data(reds, ['xx'], gain_scatter=.0099999)
        w = dict([(k,1.) for k in d.keys()])
        sol0 = dict([(k,np.ones_like(v)) for k,v in gains.items()])
        sol0.update(info.compute_ubls(d,sol0))
        #sol0 = info.logcal(d)
        #for k in sol0: sol0[k] += .01*capo.oqe.noise(sol0[k].shape)
        meta, sol = info.lincal(d, sol0)
        for i in xrange(NANTS):
            self.assertEqual(sol[(i,'x')].shape, (10,10))
        for bls in reds:
            ubl = sol[bls[0]+('xx',)]
            self.assertEqual(ubl.shape, (10,10))
            for bl in bls:
                d_bl = d[bl+('xx',)]
                mdl = sol[(bl[0],'x')] * sol[(bl[1],'x')].conj() * ubl
                np.testing.assert_almost_equal(np.abs(d_bl), np.abs(mdl), 10)
                np.testing.assert_almost_equal(np.angle(d_bl*mdl.conj()), 0, 10)
        
        

if __name__ == '__main__':
    unittest.main()
