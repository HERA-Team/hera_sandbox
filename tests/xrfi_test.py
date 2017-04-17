import unittest, glob
import capo
import numpy as np
import pylab as plt

SIZE = 100
VERBOSE = True
PLOT = True

FILES = {
    'paper': glob.glob('xrfi_data/paper/chisq0*.npz'),
    'hera': glob.glob('xrfi_data/hera/chisq0*.npz'),
}

def get_accuracy(f, rfi, verbose=VERBOSE):
    correctly_flagged = np.average(f[rfi])
    m = f.copy()
    m[rfi] = 0
    false_positive = float(np.sum(m)) / (m.size - len(rfi[0]))
    if verbose:
        print '\t Found RFI: %1.3f\n\t False Positive: %1.3f' % (correctly_flagged, false_positive)
    return correctly_flagged, false_positive

def plot_waterfall(data, f, mx=10, drng=10, mode='lin'):
    if not PLOT: return
    plt.subplot(121); capo.plot.waterfall(data, mode='lin', mx=10, drng=10); plt.colorbar()
    plt.subplot(122); capo.plot.waterfall(f, mode='lin', mx=10, drng=10); plt.colorbar()
    plt.show()

def plot_result(f, rfi):
    if not PLOT: return
    plt.plot(rfi[0], rfi[1], 'ko')
    fi = np.where(f)
    plt.plot(fi[0], fi[1], 'r.')
    plt.show()

class Template(unittest.TestCase):
    rfi_gen = None # Need to override this for each TestCase, usually in setUp
    def _run_test(self, func, correct_flag, false_positive, nsig=4):
        for data,rfi in self.rfi_gen():
            f = func(data)
            if VERBOSE: print self.__class__, func.__name__
            plot_waterfall(data, f)
            f = np.where(f > nsig, 1, 0)
            cf,fp = get_accuracy(f, rfi)
            self.assertGreater(cf, correct_flag)
            self.assertLess(fp, false_positive)
    def test_detrend_deriv(self):
        self._run_test(capo.xrfi.detrend_deriv, .9, .1, nsig=4)
    def test_detrend_medfilt(self):
        self._run_test(capo.xrfi.detrend_medfilt, .99, .01, nsig=4)
    def test_detrend_medminfilt(self):
        self._run_test(capo.xrfi.detrend_medminfilt, .97, .05, nsig=6)
    def test_xrfi_simple(self):
        self._run_test(capo.xrfi.xrfi_simple, .99, .10, nsig=.5)
    def test_xrfi(self):
        self._run_test(capo.xrfi.xrfi, .99, .01, nsig=.5)

class TestSparseScatter(Template):
    def setUp(self):
        RFI = 50
        NTRIALS = 10
        NSIG = 10
        def rfi_gen():
            for i in xrange(NTRIALS):
                data = np.random.normal(size=(SIZE,SIZE))
                rfi = (np.random.randint(SIZE, size=RFI), np.random.randint(SIZE, size=RFI))
                data[rfi] = NSIG
                yield data, rfi
            return
        self.rfi_gen = rfi_gen

class TestDenseScatter(Template):
    def setUp(self):
        RFI = 1000
        NTRIALS = 10
        NSIG = 10
        def rfi_gen():
            for i in xrange(NTRIALS):
                data = np.random.normal(size=(SIZE,SIZE))
                rfi = (np.random.randint(SIZE, size=RFI), np.random.randint(SIZE, size=RFI))
                data[rfi] = NSIG
                yield data, rfi
            return
        self.rfi_gen = rfi_gen

class TestCluster(Template):
    def setUp(self):
        RFI = 10
        NTRIALS = 10
        NSIG = 10
        def rfi_gen():
            for i in xrange(NTRIALS):
                data = np.random.normal(size=(SIZE,SIZE))
                x,y = (np.random.randint(SIZE-1, size=RFI), np.random.randint(SIZE-1, size=RFI))
                x = np.concatenate([x,x,x+1,x+1])
                y = np.concatenate([y,y+1,y,y+1])
                rfi = (np.array(x), np.array(y))
                data[rfi] = NSIG
                yield data, rfi
            return
        self.rfi_gen = rfi_gen
    def test_xrfi_simple(self):
        self._run_test(capo.xrfi.xrfi_simple, .39, .10, nsig=.5) # yuck

class TestLines(Template):
    def setUp(self):
        RFI = 3
        NTRIALS = 10
        NSIG = 10
        def rfi_gen():
            for i in xrange(NTRIALS):
                data = np.random.normal(size=(SIZE,SIZE))
                x,y = (np.random.randint(SIZE, size=RFI), np.random.randint(SIZE, size=RFI))
                mask = np.zeros_like(data)
                mask[x] = 1
                mask[:,y] = 1
                data += mask * NSIG
                yield data, np.where(mask)
            return
        self.rfi_gen = rfi_gen
    def test_detrend_deriv(self):
        self._run_test(capo.xrfi.detrend_deriv, .0, .10, nsig=4) # awful
    def test_xrfi_simple(self):
        self._run_test(capo.xrfi.xrfi_simple, .75, .10, nsig=.5) # not great
    def test_xrfi(self):
        self._run_test(capo.xrfi.xrfi, .98, .01, nsig=.5) # not great

class TestBackground(Template):
    def setUp(self):
        RFI = 50
        NTRIALS = 10
        NSIG = 10
        def rfi_gen():
            for i in xrange(NTRIALS):
                sin_t = np.sin(np.linspace(0,2*np.pi,SIZE)); sin_t.shape = (-1,1)
                sin_f = np.sin(np.linspace(0,4*np.pi,SIZE)); sin_t.shape = (1,-1)
                data = sin_t * sin_f + np.random.normal(size=(SIZE,SIZE))
                rfi = (np.random.randint(SIZE, size=RFI), np.random.randint(SIZE, size=RFI))
                data[rfi] = NSIG
                yield data, rfi
            return
        self.rfi_gen = rfi_gen

class TestHERA(Template):
    def setUp(self):
        def rfi_gen():
            for f in FILES['hera']:
                data = np.load(f)['chisq']
                rfi = np.where(capo.xrfi.xrfi(data)) # XXX actual answers?
                yield data, rfi
            return
        self.rfi_gen = rfi_gen
            
class TestPAPER(Template):
    def setUp(self):
        def rfi_gen():
            for f in FILES['paper']:
                data = np.load(f)['chisq']
                rfi = np.where(capo.xrfi.xrfi(data)) # XXX actual answers?
                yield data, rfi
            return
        self.rfi_gen = rfi_gen
            

# TODO: noise tilts
# TODO: faint RFI
# TODO: combination of everything


if __name__ == '__main__':
    unittest.main()
