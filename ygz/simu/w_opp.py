#! /usr/bin/env python
import aipy as a, numpy as np, capo as C, pylab as plt
from scipy import signal
PLOT = True
#@p.ion()
#fqs = np.linspace(.1,.2,203)
class OppSolver:
    '''uses convolution to compute Opp for two visibilities, thus need to obtain two
    time series'''
    def __init__(self, fq=.15, dT=0.01, cal='psa6622_v003', beam='PAPER'):
        print 'Initializing Oppsolver'
        self.fq = fq
        self.cal = cal
        self.aa = a.cal.get_aa(self.cal, np.array([fq]))
        self.REDNORM = 1.
        self.h = a.healpix.HealpixMap(nside=64)
        tx,ty,tz = self.h.px2crd(np.arange(self.h.map.size), ncrd=3)
        #Create equatorial coordinates of the first frame T0
        self.top0 = np.array([tx,ty,tz], dtype=tx.dtype)
        self.beam = beam
        self.dT = dT
        self.k = -2j*np.pi*self.fq #to multiply fringes
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
        
    def prepare_coord(self):
        """prepares bms in memory, need more than 5MB * T1.size memory"""
        self.T0 = 2456681.50+self.dT
        self.T1 = np.arange(2456681.4,2456681.6+self.dT, self.dT)
        self.aa.set_jultime(self.T0)
        m = np.linalg.inv(self.aa.eq2top_m)
        ex,ey,ez = np.dot(m, self.top0)
        self.eq = np.array([ex,ey,ez], dtype=ex.dtype)
        self.bms = np.zeros((self.T1.size, ex.size), dtype=np.complex)
        self.tops = np.zeros((self.T1.size, 3, ex.size))
        for i, t1 in enumerate(self.T1):
            #print t1
            self.aa.set_jultime(t1)
            m = self.aa.eq2top_m
            tx,ty,tz = np.dot(m, self.eq)
            if self.beam == 'HERA':
                bm = self.get_hera_beam((tx,ty,tz), pol='I')**2
            elif self.beam == 'PAPER':
                bm = self.aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/np.abs(tz)#*np.abs(tzsave)
            #bm = np.ones_like(tx)
            #bm = np.where(tz > 0, bm, 0)
            bm = np.where(tz > 0.001, bm, 0)
            self.bms[i] = bm
            self.tops[i] = np.asarray([tx,ty,tz])
        self.tops = self.tops.transpose((1,0,2))  #put tops in (3, time, space)


        self.REDNORM,self.Tac_err = self.w_opp((103,26),(103,26))
        print 'self.REDNORM, self.Tac_err= ', self.REDNORM, self.Tac_err

    def w_opp(self, bl1=None,bl2=None, bl1coords=None, bl2coords=None):
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
        fng1 = np.exp(self.k*bl1_prj)
        fng2 = np.exp(self.k*bl2_prj)
        V2 = self.bms * fng2
        V1 = self.bms * fng1

        # datashape = (self.T1.size, self.tops[0][0].size)
        # V1 = np.zeros(datashape, dtype=np.complex)
        # V2 = np.zeros(datashape, dtype=np.complex)
        # for i,t1 in enumerate(self.T1):
        #     tx,ty,tz = self.tops[i]
        #     bm = self.bms[i]
        #     bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
        #     bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
        #     fng1 = np.exp(-2j*np.pi*bl1_prj*self.fq)
        #     fng2 = np.exp(-2j*np.pi*bl2_prj*self.fq)
        #     bm_fng2 = bm * fng2
        #     bm_fng1 = bm * fng1
        #     V1[i] = bm_fng1
        #     V2[i] = bm_fng2
        #V1,V2 = np.array(V1), np.array(V2) 
        _V1,_V2 = np.fft.fft(V1,axis=0),np.fft.fft(V2,axis=0)
        #import IPython; IPythonp.embed()
        res = np.fft.ifftshift(np.fft.ifft(np.sum(_V2*np.conj(_V1),axis=1)))
        #res = np.fft.fftshift(np.sum(_V2*np.conj(_V1),axis=1))
        ###################
        res = res/self.REDNORM
        ###################
        #import IPython; IPythonp.embed()

        maxind = np.argmax(np.abs(res))
        maxres = res[maxind]
        T1ac = -self.T0+self.T1[maxind]
        # print '############## OPP RESULT for', bl1, bl2, '#####################'
        # print 'max, abs(max), dT(max)'
        # print maxres,maxres, T1ac
        return maxres,T1ac