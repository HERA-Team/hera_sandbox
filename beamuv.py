import aipy as a, numpy as n, pylab as p

def coeffs_from_file(infile):
    _coeffs = n.load(infile)
    return _coeffs['coeffs']

class BeamUV(a.phs.Beam):
    '''
    A primary beam model which is stored as a set of CLEAN components of the Fourier Transform of the Beam.  Once Fourier transformed, the coordinates are l,m direction cosines on the sky.
    '''
    def __init__(self,_coeffs,freqs,pol='y',size=64,res=0.4):
        self.freqs = freqs
        self.size = size
        self.res = res
        center = size/2
        offset = _coeffs.shape[0]/2
        _beam = n.zeros((size,size))
        print _beam[center-offset:center+offset+1,center-offset:center+offset+1].shape
        _beam[center-offset:center+offset+1,center-offset:center+offset+1] = _coeffs
        beam = n.fft.ifft2(_beam)*(_beam.size)
        if pol == 'x':
            beam = n.rot90(beam)
        self.im = n.sqrt(n.abs(beam))

    def xyz2lm(self,c1,c2,c3):
        '''
        Convert topocentric coordinates to direction cosines.
        '''
        l,m = c1,c2
        return l,m
    def altaz2lm(self,c1,c2):
        '''
        Convert alt/az coordinates to l,m.  Input is alt, az.
        '''
        x,y,z = a.coord.azalt2top((c2,c1))
        l,m = x,y
        return l,m
    def response(self,c1,c2,c3=None):
        '''
        Calculate the beam response in a given direction.  If two coordinates ar given, they are treated as equatorial; if three, topcentric.
        '''
        if c3 == None:
            print "Equatorial support not yet enabled."
        else:
            l,m = self.xyz2lm(c1,c2,c3)
            #l,m = (self.size*l*self.res)+(self.size/2),(self.size*m*self.res)+(self.size/2)
            l,m = (self.size*l*self.res),(self.size*m*self.res)
            #You need better interpolation.
            cl,fl = n.ceil(l).astype(int), n.floor(l).astype(int)
            cm,fm = n.ceil(m).astype(int), n.floor(m).astype(int)
            wl,wm = l%1,m%1
            r00 = self.im[fl,fm]
            r01 = self.im[fl,cm]
            r10 = self.im[cl,fm]
            r11 = self.im[cl,cm]
            response = (r00*(1-wl)*(1-wm) + r01*(1-wl)*wm + r10*wl*(1-wm) + r11*wl*wm)
            #print l,m,wl,wm,response**2
            return (n.ones_like(self.freqs) * response)
    def _plot_horizon_(self):
        h_az = n.arange(0,2*n.pi+0.5,.05)
        h_alt = n.zeros_like(h_az)
        h_l,h_m = self.altaz2lm(h_alt,h_az)
        h_l,h_m = (self.size*h_l*self.res)+(self.size/2),(self.size*h_m*self.res)+(self.size/2)
        p.plot(h_l,h_m,color='k')
    def _show_(self,power=False):
        im = a.img.recenter(self.im,(-self.size/2,-self.size/2))
        if power==True:
            p.imshow(im**2,extent=[0,self.size,0,self.size],interpolation='nearest')
        else:
            p.imshow(im,extent=[0,self.size,0,self.size],interpolation='nearest')
        self._plot_horizon_()
        p.colorbar(shrink=0.5)
    def showtrack(self,c1,c2,c3=None):
        if c3 == None:
            print "Equatorial support not yet enabled."
        else:
            l,m = self.xyz2lm(c1,c2,c3)
            l,m = (self.size*l*self.res)+(self.size/2),(self.size*m*self.res)+(self.size/2)
        self._show_(power=True)
        self._plot_horizon_()
        p.plot(l,m,'.',color='k')
