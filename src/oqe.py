import numpy as np, aipy, random

DELAY = False

def noise(size):
    sig = 1./np.sqrt(2)
    return np.random.normal(scale=sig, size=size) + 1j*np.random.normal(scale=sig, size=size)

def cov(d1, w1, d2=None, w2=None):
    if d2 is None: d2,w2 = d1.conj(),w1
    d1sum,d1wgt = (w1*d1).sum(axis=1), w1.sum(axis=1)
    d2sum,d2wgt = (w2*d2).sum(axis=1), w2.sum(axis=1)
    x1,x2 = d1sum / np.where(d1wgt > 0,d1wgt,1), d2sum / np.where(d2wgt > 0,d2wgt,1)
    x1.shape = (-1,1)
    x2.shape = (-1,1)
    d1x = d1 - x1
    d2x = d2 - x2
    C = np.dot(w1*d1x,(w2*d2x).T)
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
    def __init__(self, dsets=None, wgts=None, lsts=None, conj=None, npzfile=None, lmin=0):
        self.x, self.w = {}, {}
        self.clear_cache()
        self.lmin = lmin
        if not npzfile is None: self.from_npz(npzfile)
        elif not dsets is None: self.set_data(dsets, wgts=wgts, conj=conj)
    def flatten_data(self, data):
        if data is None: return None
        d = {}
        for k in data:
            for bl in data[k]:
                for pol in data[k][bl]:
                    key = (k, bl, pol)
                    d[key] = data[k][bl][pol]
        return d
    def set_data(self, dsets, wgts=None, conj=None):
        if type(dsets.values()[0]) == dict:
            dsets,wgts = self.flatten_data(dsets), self.flatten_data(wgts)
        self.x, self.w = {}, {}
        for k in dsets:
            self.x[k] = dsets[k].T
            try: self.w[k] = wgts[k].T
            except(TypeError): self.w[k] = np.ones_like(self.x[k])
            try:
                if conj[k[1]]: self.x[k] = np.conj(self.x[k])
            except(TypeError,KeyError): pass
    def add_data(self, dsets, wgts=None, conj=None):
        if type(dsets.values()[0]) == dict:
            dsets,wgts = self.flatten_data(dsets), self.flatten_data(wgts)
        for k in dsets:
            self.x[k] = dsets[k].T
            try: self.w[k] = wgts[k].T
            except(TypeError): self.w[k] = np.ones_like(self.x[k])
            try:
                if conj[k[1]]: self.x[k] = np.conj(self.x[k])
            except(TypeError,KeyError): pass
    def lst_align(self, lsts, dsets, wgts=None):
        for k in lsts: #orders LSTs
            order = np.argsort(lsts[k])
            lsts[k] = lsts[k][order]
        numkeys = len(lsts.keys())
        i=0 
        while i < numkeys-1: #aligns LSTs
            if i==0: lsts_final = np.intersect1d(lsts[lsts.keys()[i]],lsts[lsts.keys()[i+1]]) #XXX LSTs much match exactly
            else: lsts_final = np.intersect1d(lsts_final,lsts[lsts.keys()[i+1]])
            i += 1
        if numkeys == 1: lsts_final = lsts[lsts.keys()[0]]
        ind = {}
        for k in lsts:
            ind[k] = lsts[k].searchsorted(lsts_final)
        for k in dsets:
            dsets[k] = dsets[k][ind[k[0]]]
            if wgts: wgts[k] = wgts[k][ind[k[0]]]
        return lsts[k[0]][ind[k[0]]], dsets, wgts #lsts computed from last k but it doesn't matter

    def clear_cache(self, keys=None):
        if keys is None: self._C, self._Ctrue, self._iC = {}, {}, {}
        else:
            for k in keys:
                try: del(self._C[k])
                except(KeyError): pass
                try: del(self._Ctrue[k])
                except(KeyError): pass
                try: del(self._iC[k])
                except(KeyError): pass
    def C(self, k):
        if not self._C.has_key(k): # defaults to true covariance matrix
            self.set_C({k:cov(self.x[k], self.w[k])})
            self._Ctrue[k] = self._C[k] # save computing this later, if we are making it now
        return self._C[k]
    def set_C(self, d):
        self.clear_cache(d.keys())
        for k in d: self._C[k] = d[k]
    def Ctrue(self, k):
        if not self._Ctrue.has_key(k): self._Ctrue[k] = cov(self.x[k], self.w[k]) # must be actual covariance, no overwriting
        return self._Ctrue[k]
    def iC(self, k):
        if not self._iC.has_key(k):
            C = self.C(k)
            U,S,V = np.linalg.svd(C.conj()) # conj in advance of next step
            S += self.lmin # ensure invertibility
            self.set_iC({k:np.einsum('ij,j,jk', V.T, 1./S, U.T)})
        return self._iC[k]
    def set_iC(self, d):
        for k in d: self._iC[k] = d[k]
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
    def gen_bl_boots(self, nboots, ngps=5):
        _bls = {}
        for k in self.x: _bls[k[1]] = None
        for boot in xrange(nboots):
            bls = _bls.keys()[:]
            random.shuffle(bls)
            gps = [bls[i::ngps] for i in range(ngps)]
            gps = [[random.choice(gp) for bl in gp] for gp in gps]
            yield gps
        return
    def q_hat(self, k1, k2, use_cov=True, use_fft=True):
        nchan = self.x[k1].shape[0]
        if use_cov: iC1,iC2 = self.iC(k1), self.iC(k2)
        else:
            #iC1 = np.linalg.inv(self.C(k1) * np.identity(nchan))
            #iC2 = np.linalg.inv(self.C(k2) * np.identity(nchan))
            iC1 = iC2 = np.identity(nchan)
        if use_fft:
            iC1x, iC2x = np.dot(iC1, self.x[k1]), np.dot(iC2, self.x[k2])
            _iC1x, _iC2x = np.fft.fft(iC1x.conj(), axis=0), np.fft.fft(iC2x.conj(), axis=0)
            return np.fft.fftshift(_iC1x,axes=0).conj() * np.fft.fftshift(_iC2x,axes=0)
        else: # slow, used to explicitly cross-check fft code
            q = []
            for i in xrange(nchan):
                Q = get_Q(i,nchan)
                iCQiC = np.einsum('ab,bc,cd', iC1.T.conj(), Q, iC2) # C^-1 Q C^-1
                qi = np.sum(self.x[k1].conj() * np.dot(iCQiC,self.x[k2]), axis=0)
                q.append(qi)
            return np.array(q)
    def get_F(self, k1, k2, use_cov=True):
        nchan = self.x[k1].shape[0]
        F = np.zeros((nchan,nchan), dtype=np.complex)
        if use_cov:
            iC1,iC2 = self.iC(k1), self.iC(k2)
            Ctrue1, Ctrue2 = self.Ctrue(k1), self.Ctrue(k2)
        else:
            #iC1 = np.linalg.inv(self.C(k1) * np.identity(nchan))
            #iC2 = np.linalg.inv(self.C(k2) * np.identity(nchan))
            iC1 = iC2 = np.identity(nchan, dtype=F.dtype)
            Ctrue1 = Ctrue2 = np.identity(nchan, dtype=F.dtype) # XXX why do this
        #Ctrue1, Ctrue2 = self.Ctrue(k1), self.Ctrue(k2)
        if False: # This is for the "true" Fisher matrix
            CE1, CE2 = {}, {}
            for ch in xrange(nchan):
                Q = get_Q(ch,nchan)
                CE1[ch] = np.dot(Ctrue1, np.dot(iC1, np.dot(Q, iC2))) # C1 Cbar1^-1 Q Cbar2^-1
                CE2[ch] = np.dot(Ctrue2, np.dot(iC2, np.dot(Q, iC1))) # C2 Cbar2^-1 Q Cbar1^-1
                #CE1[ch] = np.einsum('ab,bc,cd,de', self.Ctrue(k1), iC1, Q, iC2) # slow
                #CE2[ch] = np.einsum('ab,bc,cd,de', self.Ctrue(k2), iC2, Q, iC1) # slow
            #import IPython; IPython.embed()
            for i in xrange(nchan):
                for j in xrange(nchan):
                    F[i,j] += np.einsum('ij,ji', CE1[i], CE2[j]) # C E C E
        else: # This is for the "effective" matrix s.t. W=MF and p=Mq
            iCQ1,iCQ2 = {}, {}
            for ch in xrange(nchan): # this loop is nchan^3
                Q = get_Q(ch,nchan)
                iCQ1[ch] = np.dot(iC1,Q) #C^-1 Q
                iCQ2[ch] = np.dot(iC2,Q) #C^-1 Q
            for i in xrange(nchan): # this loop goes as nchan^4
                for j in xrange(nchan):
                    F[i,j] += np.einsum('ij,ji', iCQ1[i], iCQ2[j]) #C^-1 Q C^-1 Q 
        return F
    def get_MW(self, F, mode='F^-1'):
        modes = ['F^-1', 'F^-1/2', 'I', 'L^-1']; assert(mode in modes)
        if mode == 'F^-1':
            U,S,V = np.linalg.svd(F)
            M = np.einsum('ij,j,jk', V.T, 1./S, U.T)
        elif mode == 'F^-1/2':
            U,S,V = np.linalg.svd(F)
            M = np.einsum('ij,j,jk', V.T, 1./np.sqrt(S), U.T)
        elif mode == 'I':
            M = np.identity(F.shape[0], dtype=F.dtype)
        else:
            #Cholesky decomposition to get M
            order = np.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1]) # XXX needs generalizing
            iorder = np.argsort(order)
            F_o = np.take(np.take(F,order, axis=0), order, axis=1)
            L_o = np.linalg.cholesky(F_o)
            U,S,V = np.linalg.svd(L_o.conj())
            M_o = np.dot(np.transpose(V), np.dot(np.diag(1./S), np.transpose(U)))
            M = np.take(np.take(M_o,iorder, axis=0), iorder, axis=1)
        W = np.dot(M, F)
        norm  = W.sum(axis=-1); norm.shape += (1,)
        M /= norm; W = np.dot(M, F)
        return M,W
    def p_hat(self, M, q, scalar=1.):
        return np.dot(M, q) * scalar

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
