#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as plt
from scipy import signal
PLOT = True
#@p.ion()
#fqs = n.linspace(.1,.2,203)
class OppSolver:
    '''uses convolution to compute Opp for two visibilities, thus need to obtain two
    time series'''
    def __init__(self, fq=.15, cal='psa6240_v003'):
        print 'Initializing Oppsolver'
        self.fq = fq
        self.cal = cal
        self.aa = a.cal.get_aa(self.cal, n.array([fq]))
        self.REDNORM = 1.
        self.h = a.healpix.HealpixMap(nside=64)
        tx,ty,tz = self.h.px2crd(n.arange(self.h.map.size), ncrd=3)
        #Create equatorial coordinates of the first frame T0
        self.top0 = n.array([tx,ty,tz], dtype=tx.dtype)
        self.prepare_coord()
        
    def prepare_coord(self):
        """prepares bms in memory, need more than 5MB * T1.size memory"""
        self.T0 = 2456681.501
        self.T1 = n.arange(2456681.3,2456681.701,0.001)
        self.aa.set_jultime(self.T0)
        m = n.linalg.inv(self.aa.eq2top_m)
        ex,ey,ez = n.dot(m, self.top0)
        self.eq = n.array([ex,ey,ez], dtype=ex.dtype)
        self.bms = []
        self.tops = []
        for t1 in self.T1:
            print t1
            self.aa.set_jultime(t1)
            m = self.aa.eq2top_m
            tx,ty,tz = n.dot(m, self.eq)
            bm = self.aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/n.abs(tz)#*n.abs(tzsave)
            #bm = n.ones_like(tx)
            #bm = n.where(tz > 0, bm, 0)
            bm = n.where(tz > 0.001, bm, 0)
            self.bms.append(bm)
            self.tops.append((tx,ty,tz))
        self.REDNORM,self.Tac_err = self.w_opp((103,26),(103,26))
        print 'self.REDNORM, self.Tac_err= ', self.REDNORM, self.Tac_err

    def w_opp(self, bl1,bl2):
        #h = a.healpix.HealpixMap(nside=64
        bl1x, bl1y, bl1z = self.aa.get_baseline(bl1[0],bl1[1],'z')
        bl2x, bl2y, bl2z = self.aa.get_baseline(bl2[0],bl2[1],'z')
        V1,V2 = [],[]
        #print 'Computing Opp'
        for i,t1 in enumerate(self.T1):
            tx,ty,tz = self.tops[i]
            bm = self.bms[i]
            bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
            bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
            fng1 = n.exp(-2j*n.pi*bl1_prj*self.fq)
            fng2 = n.exp(-2j*n.pi*bl2_prj*self.fq)
            bm_fng2 = bm * fng2
            bm_fng1 = bm * fng1
            V1.append(bm_fng1)
            V2.append(bm_fng2)
        V1,V2 = n.array(V1), n.array(V2) 
        _V1,_V2 = n.fft.fft(V1,axis=0),n.fft.fft(V2,axis=0)
        #import IPython; IPython.embed()
        res = n.fft.ifftshift(n.fft.ifft(n.sum(_V2*n.conj(_V1),axis=1)))
        #res = n.fft.fftshift(n.sum(_V2*n.conj(_V1),axis=1))
        ###################
        res = res/self.REDNORM
        ###################
        #import IPython; IPython.embed()

        maxind = n.argmax(n.abs(res))
        maxres = res[maxind]
        T1ac = -self.T0+self.T1[maxind]
        # print '############## OPP RESULT for', bl1, bl2, '#####################'
        # print 'max, abs(max), dT(max)'
        # print maxres,maxres, T1ac
        return maxres,T1ac