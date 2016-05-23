import numpy as np, aipy

DELAY = False

def cov(d1, w1, d2=None, w2=None):
    if d2 is None: d2,w2 = d1.conj(),w1
    d1sum,d1wgt = d1.sum(axis=1), w1.sum(axis=1)
    d2sum,d2wgt = d2.sum(axis=1), w2.sum(axis=1)
    x1,x2 = d1sum / np.where(d1wgt > 0,d1wgt,1), d2sum / np.where(d2wgt > 0,d2wgt,1)
    x1.shape = (-1,1)
    x2.shape = (-1,1)
    d1x = d1 - x1
    d2x = d2 - x2
    C = np.dot(d1x,d2x.T)
    W = np.dot(w1,w2.T)
    return C / np.where(W > 1, W-1, 1)

def get_Q(mode, n_k, window='none'): #encodes the fourier transform from freq to delay
    if not DELAY:
        _m = np.zeros((n_k,), dtype=np.complex)
        _m[mode] = 1. #delta function at specific delay mode
        m = np.fft.fft(np.fft.ifftshift(_m)) * aipy.dsp.gen_window(n_k, window) #FFT it to go to freq
        Q = np.einsum('i,j', m, m.conj()) #dot it with its conjugate
        return Q
    else:
        # XXX need to have this depend on window
        Q = np.zeros_like(C)
        Q[mode,mode] = 1
        return Q

class DataSet:
    def __init__(self, dsets=None, wgts=None, lsts=None, conj=None, npzfile=None):
        if not npzfile is None: self.from_npz(npzfile)
        else:
            self.x, self.w = {}, {}
            self.dsets = dsets.keys()
            self.bls = dsets.values()[0].keys()
            self.pols = dsets.values()[0].values()[0].keys()
            self.lsts = lsts
            for k in self.dsets:
                for bl in self.bls:
                    for pol in self.pols:
                        key = (k, bl, pol)
                        self.x[key] = dsets[k][bl][pol].T
                        self.w[key] = wgts[k][bl][pol].T
                        if conj[bl]: self.x[key] = np.conj(self.x[key])
            #self.gen_covs()
            #self.gen_icovs()
    def lst_align(self, lst1, lst2, lstres=.001):
        # XXX super ugly brute force
        i1,i2 = [], []
        for i,L1 in enumerate(lst1):
            for j,L2 in enumerate(lst2):
                match = (np.abs(L1-L2) < lstres)
                if match: break
            if match:
                i1.append(i)
                i2.append(j)
        return np.array(i1), np.array(i2)
    def gen_covs(self):
        self.C = {}
        for k in self.x: self.C[k] = cov(self.x[k], self.w[k])
    def gen_icovs(self, lmin=0):
        self.iC = {}
        for k in self.C:
            C = self.C[k]
            U,S,V = np.linalg.svd(C.conj()) # conj in advance of next step
            S += lmin # ensure invertibility
            self.iC[k] = np.einsum('ij,j,jk', V.T, 1./S, U.T)
    def to_npz(self, filename):
        data = {}
        for k in self.x:
            sk = str(k)
            data['x '+sk] = self.x[k]
            try: data['C '+sk] = self.C[k]
            except(AttributeError): pass
            try: data['iC '+sk] = self.iC[k]
            except(AttributeError): pass
        for k in self.lsts:
            data['lst '+k] = self.lsts[k]
        np.savez(filename, **data)
    def from_npz(self, filename):
        npz = np.load(filename)
        self.x = {}
        for k in [f for f in npz.files if f.startswith('x')]: self.x[eval(k[2:])] = npz[k]        
        self.C = {}
        for k in [f for f in npz.files if f.startswith('C')]: self.C[eval(k[2:])] = npz[k]        
        self.iC = {}
        for k in [f for f in npz.files if f.startswith('iC')]: self.iC[eval(k[3:])] = npz[k]        
        self.lsts = {}
        for k in [f for f in npz.files if f.startswith('lst')]: self.lsts[k[4:]] = npz[k]        
        dsets, bls, pols = {}, {}, {}
        for k in self.x:
            dsets[k[0]] = None
            bls[k[1]] = None
            pols[k[2]] = None
        self.dsets = dsets.keys()
        self.bls = bls.keys()
        self.pols = pols.keys()

'''
def oqe(dsets, conj, chans, ):
    nchan = chans.size
    bls_master = dsets.values()[0].keys()
    nbls = len(bls_master)
    Q = [get_Q(i,nchan) for i in xrange(nchan)]
    
    I,_I,_Ix = {},{},{}
    C,_C,_Cx = {},{},{}
    for k in dsets:
        I[k],_I[k],_Ix[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in x[k]:
            C[k][bl] = cov(x[k][bl])
            #C[k][bl] = covs[k][str(bl+(POL,))][120:141,120:141]
            #if conj[bl]: C[k][bl] = C[k][bl].conj()
            I[k][bl] = n.identity(C[k][bl].shape[0])
            U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
            _Ix[k][bl] = x[k][bl].copy()
    for boot in xrange(opts.nboot):
        print '%d / %d' % (boot+1,opts.nboot)
        bls = bls_master[:]
        if True: #shuffle and group baselines for bootstrapping
            if not SAMPLE_WITH_REPLACEMENT:
                random.shuffle(bls)
                bls = bls[:-5] # XXX
            else: #sample with replacement
                bls = [random.choice(bls) for bl in bls]
            gps = [bls[i::NGPS] for i in range(NGPS)]
            gps = [[random.choice(gp) for bl in gp] for gp in gps]
        else: #assign each baseline its own group
            gps = [bls[i::NGPS] for i in range(NGPS)]
        bls = [bl for gp in gps for bl in gp]
        #print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])    
        _Iz,_Isum,_IsumQ = {},{},{}
        _Cz,_Csum,_CsumQ = {},{},{}
        print "   Getting C"
        for k in days:
            _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
            _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
            for i,gp in enumerate(gps): #sum things up over the groups
                _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
                _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
                _Isum[k][i] = sum([_I[k][bl] for bl in gp])
                _Csum[k][i] = sum([_C[k][bl] for bl in gp])
                _IsumQ[k][i] = {}
                _CsumQ[k][i] = {}
                if DELAY: #this is much faster
                    _Iz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                    _Cz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Cz[k][i], axis=0), axes=0)
                    #XXX need to take fft of _Csum, _Isum here
                for ch in xrange(nchan): #XXX this loop makes computation go as nchan^3
                    _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                    _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch]) #C^-1 Q
'''
