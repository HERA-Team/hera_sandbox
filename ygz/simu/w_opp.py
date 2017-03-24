#! /usr/bin/env python
import aipy as a, numpy as np, capo as C, pylab as plt
from scipy import signal 
from scipy.interpolate import interp1d
from skimage.restoration import unwrap_phase
PLOT = True
#@p.ion()
#fqs = np.linspace(.1,.2,203)
class OppSolver:
    '''uses convolution to compute Opp for two visibilities, thus need to obtain two
    time series'''
    #@profile
    def __init__(self, fqs=np.linspace(.14,.16,20), T1=np.arange(2456681.4, 2456681.6, 0.001), cal='psa6622_v003', beam='PAPER', bandpass=None):
        print 'Initializing Oppsolver'
        self.fqs = fqs
        self.cal = cal
        self.aa = a.cal.get_aa(self.cal, fqs)
        self.h = a.healpix.HealpixMap(nside=64)
        tx,ty,tz = self.h.px2crd(np.arange(self.h.map.size), ncrd=3)
        #Create equatorial coordinates of the first frame T0
        self.top0 = np.array([tx,ty,tz], dtype=tx.dtype)
        self.beam = beam
        self.T1 = T1
        self.T0 = (T1[0]+T1[-1])/2
        self.k = -2j*np.pi*self.fqs[:,np.newaxis, np.newaxis] #to multiply fringes, prepare two extra dimensions for time and space
        if bandpass is not None:
            tfqs = np.linspace(.1,.2,bandpass.size)
            tf = interp1d(tfqs, bandpass)
            self.bandpass = tf(self.fqs)[np.newaxis, :, np.newaxis]
        else:
            self.bandpass = np.ones_like(fqs)[np.newaxis, :, np.newaxis]
        if self.beam == 'HERA':
            DATADIR = '../calfiles/'
            XFILE = DATADIR+'GX4Y2H_4900_150.hmap'
            YFILE = DATADIR+'GY4Y2H_4900_150.hmap'
            self.Xh = a.map.Map(fromfits=XFILE)
            self.Yh = a.map.Map(fromfits=YFILE)
            self.Xh.set_interpol(True)
            self.Yh.set_interpol(True)

        self.prepare_coord()

    def get_hera_beam(self, ntop, pol='I'):
        X,Y,Z = ntop
        bmI = np.sqrt((self.Xh[X,Y,Z].conj()*self.Xh[X,Y,Z] + self.Yh[X,Y,Z].conj()*self.Yh[X,Y,Z])*0.5)
        return bmI
    #@profile
    def prepare_coord(self):
        """prepares bms in memory, need more than 5MB * T1.size memory"""
        self.aa.set_jultime(self.T0)
        m = np.linalg.inv(self.aa.eq2top_m)
        ex,ey,ez = np.dot(m, self.top0)
        #zex,zey,zez = np.dot(m, np.array([0,0,1]))
        self.eq = np.array([ex,ey,ez], dtype=ex.dtype)
        self.bms = np.zeros((self.T1.size, self.fqs.size, ex.size), dtype=np.complex)
        self.tops = np.zeros((self.T1.size, 3, ex.size))
        #self.zens = np.zeros((self.T1.size, 3))
        for i, t1 in enumerate(self.T1):
            #print t1
            self.aa.set_jultime(t1)
            m = self.aa.eq2top_m
            tx,ty,tz = np.dot(m, self.eq)
            #zx,zy,zz = np.dot(m, np.array([zex,zey,zez])) #location of zenith
            if self.beam == 'HERA':
                bm = self.get_hera_beam((tx,ty,tz), pol='I')**2
            elif self.beam == 'PAPER':
                bm = self.aa[0].bm_response((tx,ty,tz),pol='I')**2#/np.abs(tz)#*np.abs(tzsave)
            bm = np.where(tz > 0.001, bm, 0)
            self.bms[i] = bm
            self.tops[i] = np.asarray([tx,ty,tz])
            #self.zens[i] = np.asarray([zx,zy,zz])
        self.bms = self.bms*self.bandpass
        self.tops = self.tops.transpose((1,0,2))[:,np.newaxis,...] 
        #self.zens = self.zens.transpose((1,0))[:,np.newaxis,:,np.newaxis] 
        #put tops in (3, 1, time, space)
        #bms should be in (freq, space)


        #self.REDNORM,self.Tac_err = self.w_opp((103,26),(103,26))
        #print 'self.REDNORM, self.Tac_err= ', self.REDNORM, self.Tac_err
    #@profile
    def opp(self, bl1=None,bl2=None, bl1coords=None, bl2coords=None, rephase=0, delay=True, return_series=False, debug=False):
        #h = a.healpix.HealpixMap(nside=64
        if bl1coords:
            #convert meters to light seconds to work with fq in GHz
            bl1x, bl1y, bl1z = np.asarray(bl1coords)/a.const.len_ns * 100. 
            bl2x, bl2y, bl2z = np.asarray(bl2coords)/a.const.len_ns * 100.

        elif bl1 is not None:
            bl1x, bl1y, bl1z = self.aa.get_baseline(bl1[0],bl1[1],'z')
            bl2x, bl2y, bl2z = self.aa.get_baseline(bl2[0],bl2[1],'z')
        else:
            raise Exception("Must supply either bl1 or bl1coords")

        tx,ty,tz = self.tops
        bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
        bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z

        # zx,zy,zz = self.zens
        # bl2_zen = zx*bl2x + zy*bl2y + zz*bl2z
        # bl1_zen = zx*bl1x + zy*bl1y + zz*bl1z

        bl1_prj -= rephase/2/np.pi
        #import IPython; IPython.embed()

        # Original Code for pedagogy, new code for memoery efficiency
        if debug:
            fng1 = np.exp(self.k*bl1_prj) #(freq, time, space)
            fng2 = np.exp(self.k*bl2_prj)  #memory efficiency
            V2 = self.bms * np.exp(self.k*bl2_prj).transpose((1,0,2)) #(time, freq, space)
            V1 = self.bms * np.exp(self.k*bl1_prj).transpose((1,0,2))
            #convolve along time axis, sum over angles and frequencies, getting (time(offset), freq)
            _V1,_V2 = np.fft.fft(V1,axis=0),np.fft.fft(V2,axis=0)  
            #import IPython; IPythonp.embed()
            res = np.fft.ifftshift(np.fft.ifft(np.sum(np.mean(_V2*np.conj(_V1),axis=2), axis=1),axis=0), axes=0)
            #res = np.fft.fftshift(np.sum(_V2*np.conj(_V1),axis=1))
        else:
            if delay:
                res = np.fft.fftshift(np.fft.ifft(
                    np.sum(np.mean(
                        np.fft.fft(self.bms * np.exp(self.k*bl2_prj).transpose((1,0,2)),axis=0)*np.conj(np.fft.fft(self.bms * np.exp(self.k*bl1_prj).transpose((1,0,2)),axis=0)),axis=2), axis=1, keepdims=True),axis=0), axes=0)
        
            else: #return results as function of frequency
                res = np.fft.fftshift(np.fft.ifft(
                    np.mean(
                        np.fft.fft(self.bms * np.exp(self.k*bl2_prj).transpose((1,0,2)),axis=0)*np.conj(np.fft.fft(self.bms * np.exp(self.k*bl1_prj).transpose((1,0,2)),axis=0)),axis=2),axis=0), axes=0)

        if return_series:
            return res
        else:
            maxind = np.argmax(np.abs(res), axis=0)
            maxres = res[maxind, np.arange(res.shape[1])]
            T1ac = -self.T0+self.T1[maxind]

            return maxres,T1ac

    def opp_rephase(self, bl1=None,bl2=None, bl1coords=None, bl2coords=None, return_series=False):
        maxres, _ = self.opp(bl1=bl1,bl2=bl2, bl1coords=bl1coords, bl2coords=bl2coords, 
            rephase=0, delay=False, return_series=False)
        slope = 0
        bl1off = unwrap_phase(np.angle(maxres))
        slope, intercept = np.polyfit(self.fqs, bl1off, 1)
        return self.opp(bl1=bl1,bl2=bl2, bl1coords=bl1coords, bl2coords=bl2coords, 
            rephase=slope, delay=True, return_series=return_series)

if __name__=='__main__':
    pass