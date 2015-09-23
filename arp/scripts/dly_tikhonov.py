#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, capo as C
import sys
import omnical as O

uv = a.miriad.UV(sys.argv[-1])
t,dat,flg = C.arp.get_dict_of_uv_data(sys.argv[-1:],antstr='1_4',polstr='xx', verbose=True)
for bl in dat:
    for pol in dat[bl]:
        d,w = dat[bl][pol], n.logical_not(flg[bl][pol]).astype(n.float)
        nchan = d.shape[-1]
        alpha = 1e-1
        I = n.identity(nchan, dtype=n.complex)
        for i in xrange(d.shape[0]):
            di,wi = d[i],w[i]
            #wi[::13] = 0
            wi[::5] = 0
            di *= wi
            #compr = 203.
            compr = 511
            #dio,norm,res,m = O.calibration_omni.deconvolve_spectra2(di,wi,(compr+1)/2, correction_weight=1e-3)
            dio,norm,res,m = O.calibration_omni.deconvolve_spectra2(di,wi,(compr+1)/2, correction_weight=1e-15)
            dio2,norm,res,m = O.calibration_omni.deconvolve_spectra2(di,wi,(compr+1)/2, correction_weight=1e-6)
            dio /= di.size/compr
            dio2 /= di.size/compr
            #dio,norm,res,m = O.calibration_omni.deconvolve_spectra2(di,wi,di.size/2)
            _dio = n.fft.ifft(dio)
            _dio = n.concatenate([_dio[:compr/2],n.zeros(di.size-compr),_dio[compr/2:]])
            _di,_wi = n.fft.ifft(di), n.fft.ifft(wi)
            _dic,info = a.deconv.clean(_di,_wi,tol=1e-6); _dic += info['res']
            #di.shape = (nchan,1)
            _di.shape = (nchan,1)
            wi.shape = (1,nchan)
            #A = n.fft.ifft(n.fft.fft(n.identity(nchan,dtype=n.complex)) * wi)
            A = n.fft.fft(n.fft.ifft(n.identity(nchan,dtype=n.complex)) * wi)
            if False: p.subplot(111); C.arp.waterfall(A); p.show()
            At = A.T.conj()
            M = n.dot(At, A)
            if True:
                G = alpha * n.identity(M.shape[0], dtype=n.complex)
                G[ 0, 0] = 0
                for i in xrange(1,20):
                    G[i,i] = 0
                    G[-i,-i] = 0
                Gt = G.T.conj()
                M += n.dot(Gt,G)
                _M = n.linalg.inv(M)
            else:
                #U,S,V = n.linalg.svd(M)
                #p.plot(S); p.show()
                _M = n.linalg.pinv(M,1e-3)
            _dim = n.dot(_M, n.dot(At,_di))

            p.subplot(211)
            p.plot(n.real(di), 'k.')
            p.plot(n.real(n.fft.fft(_dic)), 'b')
            p.plot(n.real(n.fft.fft(_dim.flatten())), 'r')
            p.ylim(-.06,.06); p.xlim(0,_dic.size)
            p.subplot(212)
            p.plot(n.real(dio))
            p.plot(n.real(dio2))
            p.ylim(-.06,.06); p.xlim(0,dio2.size)
            p.show()
            if False:
                p.subplot(131); C.arp.waterfall(          M,drng=3); p.colorbar()
                p.subplot(132); C.arp.waterfall(         _M,drng=3); p.colorbar()
                p.subplot(133); C.arp.waterfall(n.dot(M,_M),drng=3); p.colorbar()
                p.show()
            if True:
                p.plot(n.real(_di ))
                p.plot(n.real(_dic))
                p.plot(n.real(_dim))
                p.plot(n.real(_dio))
                p.ylim(-.01,.01)
                p.show()
            sys.exit(0)
            
        _d = n.fft.ifft(d, axis=-1)
        p.subplot(121); C.arp.waterfall(d)
        p.subplot(122); C.arp.waterfall(_d)
        p.show()
