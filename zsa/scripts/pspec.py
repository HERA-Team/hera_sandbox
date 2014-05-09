'''Things pertaining to power spectrum analysis'''
import numpy as n, aipy as a
import capo 



class DataReader(object):
    '''
        Data reader class. Take inputs of filenames of miriad files to read, 
        aa     = AntennaArray
        antstr = Antenna string. Use aipy scripting library.
        polstr = Polarization type. Use aipy scripting library.
        chans  = Channels to choose. Use aipt scripting library.
    '''
       
    def __init__(self, filenames, aa, antstr, polstr, chans, WINDOW='blackman-harris'):
        self.filenames = filenames 
        self.aa = aa
        self.antstr = antstr
        self.polstr = polstr
        self.achans_str = chans
        self.antpos = self.aa.ant_layout
        self.WINDOW = WINDOW #name of window function.

        #These get filled in when get_uv_data is run
        #Delay specs for every baseline. Concat in time.
        self.T = {}
        #Times derived from uv files.
        self.times = []
        #lsts derived from uvfiles.
        self.lsts = []
        # _T is the array of dspecs for achans*bls X times. 
        self._T = None
        #inttime from uvfile
        self.inttime = None
        #nchans from uvfile (full bw)
        self.nchans = None
        #window used in data. 
        self.window = None
        #bls taken from uvfile.
        self.bls = None

        #wgts to use. Populated in get_wgts()
        self.wgts = None
    
        #create a bl2sep and sep2bl dictionaries.
        self.bl2sep = {}
        self.sep2bl = {}
        
        self.get_blsep()
        self.get_uv_data()
        
    def get_uv_data(self):
        '''
           Function to get dictionaries of relavant data
           By default, returns Temperature, times, lsts, delayspectra.

        '''
        if type(self.filenames)=='str': self.filenames = [self.filenames]
        
        #get parameters from uvfile.
        uv = a.miriad.UV(self.filenames[0])
        self.inttime = uv['inttime']
        self.nchans = uv['nchan']
        self.achans = a.scripting.parse_chans(self.achans_str, self.nchans)

        for filename in self.filenames:
            print 'Reading', filename
            uv = a.miriad.UV(filename)
            a.scripting.uv_selector(uv, self.antstr, self.polstr)
            freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
            afreqs = freqs.take(self.achans)
            for (crd,t,(i,j)),d,f in uv.all(raw=True):
                if len(self.times) == 0 or self.times[-1] !=t:
                    self.times.append(t)
                    self.lsts.append(uv['lst'])
                bl = a.miriad.ij2bl(i,j)
                sep = self.bl2sep[bl] #note sep is a tuple. if first index (column separation) < 0
                if sep[0] < 0:
                    d,sep = n.conj(d), -1*n.array(sep)
                        
                #convert from jansky's to temperature units. 
                d,f = d.take(self.achans), f.take(self.achans)
                self.wgts = self.get_wgts(f)
                Trms = d*capo.pspec.jy2T(afreqs)
                #take delay transform.
                self.window = a.dsp.gen_window(Trms.size, self.WINDOW)
                _Trms = n.fft.fftshift( n.fft.ifft(self.window * Trms) )
                _Wrms = n.fft.fftshift( n.fft.ifft(self.window * self.wgts) )

                self.T[bl] = self.T.get(bl, []) + [_Trms]
            
        self.bls = self.T.keys()
        for bl in self.bls:
            self.T[bl] = n.array(self.T[bl])
        self._T = n.concatenate([self.T[bl] for bl in self.bls], axis=-1).T
            
            
        #return self.T, self.times, self.lsts, self._T

    def gen_eormdl(self, nth):
        '''
            First way to generate an eor_mdl. For every nth integration make a 
            random normal noise. 
        '''
        #run through all integrations and all baselines to insert a common random
        #gaussian signal every nth integration.
        self.eormdl = {}
        for i,t in enumerate(self.times):
            for bl in self.bls: 
                if i%nth == 0:
                    eor_mdl = n.random.normal(size=len(self.achans)) * n.exp(2j*n.pi*n.random.uniform(size=len(self.achans))) * self.wgts
                    fN = n.fft.fftshift( n.fft.ifft(self.window * eor_mdl) )
                self.eormdl[bl] = self.eormdl.get(bl, []) + [fN]
        for bl in self.bls:
            self.eormdl[bl] = n.array(self.eormdl[bl])
        #Fringe rate filter noise to match the data.
        for bl in self.eormdl:
            _N = n.fft.ifft(self.eormdl[bl], axis=0)
            _N[23:] = 0
            self.eormdl[bl] = n.fft.fft(_N, axis=0)
        self.eormdl_dspec = n.concatenate([self.eormdl[bl] for bl in self.bls], axis=-1).T
            
    def gen_noise1(self,tsys=560e3,bw=100e6,nday=47,nbl=1,npol=2,tint=43): 
        '''
            Another way of generating noise using array parameters.
            nday = number of days in integration
            nbl = number of baselines types
            npl = number of polarizations
            tint = integration time (43 for compressed data)
        ''' 
        TSYS = tsys
        B = bw/self.nchans
        NDAY = nday
        NBL = nbl
        NPOL = npol
        T_INT = tint
    
        self.N1 = {}
        for i,t in enumerate(self.times):
            for bl in self.bls:
                Trms_ = n.random.normal(size=len(self.achans)) * \
                            n.exp(2j*n.pi*n.random.uniform(size=len(self.achans)))
                Trms_ *= TSYS / n.sqrt(B * T_INT * NDAY * NBL * NPOL)
                fN = n.fft.fftshift( n.fft.ifft(self.window * Trms_) )
                self.N1[bl]  = self.N1.get(bl, []) + [fN]
        for bl in self.bls:
            self.N1[bl] = n.array(self.N1[bl])
        #Fringe rate filter noise to match the data.
        for bl in self.N1:
            _N = n.fft.ifft(self.N1[bl], axis=0)
            _N[23:] = 0
            self.N1[bl] = n.fft.fft(_N, axis=0)
        self.N1_dspec = n.concatenate([self.N1[bl] for bl in self.bls], axis=-1).T
            
    def get_blsep(self):
        '''
            Creat sep2bl and bl2sep directory for baselines in antenna array.
        '''
        for ri in range(self.antpos.shape[0]):
            for ci in range(self.antpos.shape[1]):
                for rj in range(self.antpos.shape[0]):
                    for cj in range(ci, self.antpos.shape[1]):
                        if ri >= rj and ci ==cj: continue #exclude repeats
                        sep = (cj-ci, rj-ri)
                        i,j = self.antpos[ri,ci], self.antpos[rj,cj]
                        bl = a.miriad.ij2bl(i,j)
                        if i > j:
                            i,j = j,i
                            sep = (sep[0]*-1, sep[1]*-1)
                        self.bl2sep[bl] = sep
                        self.sep2bl[sep] = self.sep2bl.get(sep,[]) + [bl]

    def get_wgts(self, f, ver='default'):
        '''
            Generate weights. 
            This should be a default function for calling 
            a method of weighting. defualt is all ones.
            input f is flags (in default case).
        '''
        if ver=='default':
            return n.logical_not(f).astype(n.float)
           
        
                            









