'''module describes covariance classes used for analysis.''' 
import numpy as n
import random, collections
import aipy as a


class CoV_base(object):
    '''Covariance class. 
        input : Data matrix and bls in set. (T, bls)
        
        bls   = bls
        T     = data matrix
        nprms = number of channels/kmodes/prms
        C     = covariance of data matrix.
    '''
    def __init__(self, T, bls):
        self.bls = bls
        self.T = T
        self.nprms = T.shape[0] / len(bls)
        self.C = cov(T)
    def get_C(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        i,j = self.bls.index(bl1), self.bls.index(bl2)
        return self.C[i*self.nprms:(i+1)*self.nprms, j*self.nprms:(j+1)*self.nprms].copy()
    def get_Cinv(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        return n.linalg.inv(self.get_C(bl1,bl2))
    def get_x(self, bl):
        i = self.bls.index(bl)
        return self.T[i*self.nprms:(i+1)*self.nprms].copy()
    def get_Ck(self, k):
        inds = n.arange(len(self.bls)) * self.nprms + k
        return self.C.take(inds,axis=0).take(inds,axis=1)


class CoV(object):
    '''Covariance class. 
        input : Data matrix and bls in set. (T, bls)
        
        bls   = bls used to make data matrix.
        T     = data matrix.
        nprms = number of channels/kmodes/prms.
        C     = covariance of data matrix.
        these are inherited from COV.

        Notes:
            1. T matrix convention (nbls*chans, time)
    '''
    def __init__(self, T, bls):
        self.nbls = len(bls)
        self.bls = bls
        self.T = T
        self.origT = T
        self.nprms = T.shape[0] / self.nbls
        self.C = cov(self.T)
        self.origC = cov(self.origT)
        self._Ctot = 1
        self.SZ = self.origT.shape[0]
        self.dspecs = []
    
    def get_C(self, bl1, bl2=None, mode=False):
        '''Get covariance matrix for given baseline pair.
           Two modes :
                diag : whatever is self.C. If diagonalize has been called 
                       this C should have covariances removed. 
                orig : the original covariance matrix.
        '''
        if bl2 is None: bl2 = bl1
        i,j = self.bls.index(bl1), self.bls.index(bl2)
        if not mode : C = self.origC
        else : C = self.C
        return C[i*self.nprms:(i+1)*self.nprms, 
                        j*self.nprms:(j+1)*self.nprms].copy()
    
    def get_Cinv(self, bl1, bl2=None, mode=False):
        '''Return inverse of covariance between two bls'''
        if bl2 is None: bl2 = bl1
        return n.linalg.inv(self.get_C(bl1,bl2,mode=mode))
    
    def get_x(self, bl, mode=False):
        '''Get input data for a given baseline.
           Two modes : 
                True (1) : diag covariance applied.
                False (0): orig. no covariance applied. 
        '''
        i = self.bls.index(bl)
        if not mode : T = self.origT
        else : T = self.T
        return T[i*self.nprms:(i+1)*self.nprms].copy()
    
    def get_Ck(self, k, mode=False):
        '''Return a single k-mode for all baselines.'''
        inds = n.arange(len(self.bls)) * self.nprms + k
        if not mode : C = self.origC
        if mode : C = self.C
        return C.take(inds,axis=0).take(inds,axis=1)

    def frfilter(self, bl):
        '''Fringe rate filter the given bl.
           This should be done to the noise matrix.'''
        ift = n.fft.ifft(self.T[bl], axis=0) 
        ift[23:] = 0 # calculated by hand for fr-filter with max_fr=1, min_fr=0.
        return ift

    def normalize(self, M):
        '''Normalize a covariance matrix.'''
        d = n.diag(M); d.shape = (1,M.shape[0])
        M /= n.sqrt(d) * 2
        return M

    def apply_gain(self, C, g = 0.1):
        '''Mutiplies by the "gain" factor.
           This factor measures the linear steps taken in the 
           approximation of the cov removal.
           Note that this factor depends on the number of bls you are 
           using. 
           Creates new covariance matrix (self._C) so as to not kill the 
           original.'''
        _C = -g*C
        return _C

    def zero_diags(self, neighbors=0):
        '''Zero out the diagonals of the covariance matrices.
           zeros out diagonal for each redundant baseline pair.
           Can also zero out modes adjacent to the diagonal. 

           n :      number of adjacent diagonals to zero out.'''
        for i in xrange(neighbors+1):
            ind = n.arange(self.SZ - i)
            for b in xrange(len(self.bls_)):
                indb = ind[:-b*self.nprms]
                if neighbors==0:
                    self._C[indb, indb+b*self.nprms] = 0.
                    self._C[indb+b*self.nprms, indb] = 0.
                else:
                    self._C[indb[:-1*i], indb[i:]+b*self.nprms] = 0.
                    self._C[indb[i:]+b*self.nprms, indb[:-1*i]] = 0.
        self._C[ind, ind] = 0
        
    def diagonalize(self, gain = 0.1, neighbors=0, niters = 8):
        '''
           Estimates and removes signal covariance from diagonaliztion
           process. This is a recursive function and updates the arrays 
           self.{C,Ctot,_Ctot, _C, T}
        ''' 
        #get ready for diag..
        self.C = cov(self.T) #make covariance matrix. Will iterate over T and C
        self.C = self.normalize(self.C)
        self._C = self.apply_gain(self.C, g=gain)
        self.zero_diags(neighbors=neighbors)
        
        #self._C made in the apply gain function.
        self._C.shape = (self.nbls, self.nprms, self.nbls, self.nprms)
        sub_C = n.zeros_like(self._C)
        #choose a (i,j) baseline cross-multiple panel in the covariance matrix.
        for i in xrange(self.nbls):
            bli = self.bls_[i]
            for j in  xrange(self.nbls):
                blj = self.bls_[j]
                
                #ensure that bli and blj belong to the same group, 
                for igp in xrange(len(self.groups)):
                    if bli in self.groups[igp] and blj in self.groups[igp]:
                        gp = self.groups[igp]
                    else: gp = None
                if gp is None: continue
                
               
                #only average using bls in the same group. and then average over
                #all other panels of covariance matrix (within group) to get the
                #average signal covariance and subtract that off so that we 
                #dont get signal liss removing residual signal covariances.
                _Csum,_Cwgt = 0., 0.
                for i_ in xrange(self.nbls):
                    bli_ = self.bls_[i_]
                    if not bli_ in gp: continue#avg over bls in same gp.
                    if bli == bli_: continue#avg over other bls in same gp. 
                    for j_ in xrange(self.nbls):
                        blj_ = self.bls_[j_]
                        if not blj_ in gp: continue#avg over bls in same gp.
                        if bli_ == blj_: continue#dont avg. over pnls with bias.
                        if blj == blj_: continue#avg over othe bls in same gp.
                        
                        _Csum += self._C[i_,:,j_,:]
                        _Cwgt += 1 

                try: 
                    sub_C[i,:,j,:] = _Csum/_Cwgt
                except(ZeroDivisionError):
                    print 'Weights are zero for %d_%d'%(i_,j_)
                    sub_C[i,:,j,:] = _Csum
            
        self._C.shape = (self.nbls*self.nprms,self.nbls*self.nprms)
        sub_C.shape = (self.nbls*self.nprms,self.nbls*self.nprms)
        self._C -= sub_C
        mask = self.create_mask()
        self._C *= mask #mask out between groups
        self._C[n.arange(self.SZ), n.arange(self.SZ)] = 1 #diag =1
        self.T = n.dot(self._C,self.T)
        self._Ctot = n.dot(self._C, self._Ctot)

        if niters > 0:
            self.diagonalize(gain = 0.1, neighbors=0, niters = niters-1)
        else:
            return None

    def create_mask(self):
        '''Create mask for masking out intergroup products'''
        mask = n.zeros_like(self._C)
       
        lgp = 0 
        for i, g in enumerate(self.groups):
            lg = lgp
            lgp = len(g) * self.nprms
            mask[i*lg:(i+1)*lgp, i*lg:(i+1)*lgp] = 1.
        return mask
                    
                
    def make_gps(self, ngps=4):
        '''Makes the gps.'''
        self.bls_ = list(random.sample(self.bls, self.nbls))
        nblspg = self.nbls/ngps
        self.groups = []
        print 'Breaking %d bls into groups of %d'%(self.nbls, nblspg)
        for i in range(ngps):
            self.groups.append(self.bls_[i*nblspg:(i+1)*nblspg]) 

        leftover = self.nbls - (nblspg*ngps)
        while leftover>0:
            i = random.choice(range(ngps))
            self.groups[i].append(self.bls_[-1*leftover])
            leftover -= 1

        #Need atleast 3 unique baselines, or else get 0 divide in diag
        for i in xrange(ngps):
            self.groups[i] = random.sample(self.groups[i], 3) + [
                random.choice(self.groups[i]) for bl in self.groups[i][
                :len(self.groups[i]) - 3]]
        self.bls_ = n.sum(self.groups)

    def get_dspecs(self):
        ''' 
            Make power spectrums for all baseline pairs in bls_ for 
            the original delay spectrums.  
            Saves in self.dspecs
        '''
        for cnt,bli in enumerate(self.bls_):
            for blj in bls_[cnt:]:
                xi = get_x(blj,mode=0)
                xj = get_x(blj,mode=0)
                self.dspecs.append(n.avg(xi * xj.conj(), axis=0))

#    def get_noise_avg(self):
         
    def remove_cov_intergp(self, bli, blj, mask, in_cov=None, gain=0.1, n_k=40, kp_hor=True):
        ''' 
            Remove leakage from particular modes for a covariance pair.
        ''' 
        Ts = n.concatenate(self.get_x(bli, mode=1), self.get_x(blj, mode=1)) 
        if in_cov:
            cx = in_cov
        else:
            cx = cov(Ts)
        cx = self.normalize(cx)
        cx =  self.apply_gain(cx, gain=gain)
        mask = n.zeros_like(cx)
        n_k = self.nprms
        #delay bins inside horizon.
        if kp_hor:
            prj_ch = xrange(int(n_k/2-n_k/10),int(n_k/2+n_k/10))
            for k in prj_ch:
                mask[k] = mask[:,k] = 1
                mask[k+n_k] = mask[:,k+n_k] = 1
        ind = n.arange(n_k)
        mask[ind, ind+n_k] = mask[ind+n_k, ind] = 0
        cx *= mask
        cx[ind, ind] = cx[ind+n_k, ind+n_k] = 1
        Ts = n.dot(cx, Ts)
        cx = cov(Ts) 
         
        
       
        

    
                            
        
    
        
    

#Useful Functions
def cov(m):
    '''Because numpy.cov is stupid and casts as float.
       Returns covariance matrix (complex).
    '''
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()


def cov2(m1,m2):
    '''Because numpy.cov is stupid and casts as float.
       Returns covariance matrix between 2 data sets.
    '''
    #return n.cov(m)
    X1 = n.array(m1, ndmin=2, dtype=n.complex)
    X2 = n.array(m2, ndmin=2, dtype=n.complex)
    X1 -= X1.mean(axis=1)[(slice(None),n.newaxis)]
    X2 -= X2.mean(axis=1)[(slice(None),n.newaxis)]
    N = X1.shape[1]
    fact = float(N - 1)
    return (n.dot(X1, X2.T.conj())/fact).squeeze()


