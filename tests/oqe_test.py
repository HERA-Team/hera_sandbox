import unittest, random, glob, os
import capo
import aipy as a, numpy as np
from mpl_toolkits.basemap import Basemap
import capo.oqe as oqe

DATADIR = '/Users/aparsons/projects/eor/psa6942/omnical_v2_xtalk/lstbin000'
TESTFILE = 'test.npz'
CAL = 'psa6622_v003'
POL = 'I'
random.seed(0)

args = glob.glob(DATADIR+'/even/lst*uvGAL') + glob.glob(DATADIR+'/odd/lst*uvGAL')
dsets = {
    'even': [x for x in args if 'even' in x],
    'odd' : [x for x in args if 'odd' in x]
}
aa = a.cal.get_aa(CAL, np.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
antstr = 'cross'
lsts,data,flgs = {},{},{}
days = dsets.keys()
for k in days:
    info,data[k],flgs[k] = capo.miriad.read_files(dsets[k], antstr=antstr, polstr=POL, verbose=True)
wgts = {}
for k in flgs:
    if not wgts.has_key(k): wgts[k] = {}
    for bl in flgs[k]:
        if not wgts[k].has_key(bl): wgts[k][bl] = {}
        for pol in flgs[k][bl]:
            wgts[k][bl][pol] = np.logical_not(flgs[k][bl][pol]).astype(np.int)

class TestDataSet(unittest.TestCase):
    def test_init(self):
        ds = capo.oqe.DataSet(data,wgts,lsts,conj)
        for k in data:
            for bl in data[k]:
                for pol in data[k][bl]:
                    self.assertTrue(ds.x.has_key((k,bl,pol)))
    def test_to_from_npz(self):
        ds1 = capo.oqe.DataSet(data,wgts,lsts,conj)
        ds1.to_npz(TESTFILE)
        ds2 = capo.oqe.DataSet(npzfile=TESTFILE)
        for k in ds1.x:
            self.assertTrue(ds2.x.has_key(k))
    def tearDown(self):
        if os.path.exists(TESTFILE): os.remove(TESTFILE)
        
#ds.x
#capo.plot.waterfall(ds.x[('odd',(1,4),'I')], drng=3)
#p.show()
#ds.to_npz('test.npz')

if __name__ == '__main__':
    unittest.main()
