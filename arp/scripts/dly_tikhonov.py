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
        alpha = 1e-2
        I = n.identity(nchan, dtype=n.complex)
        for i in xrange(d.shape[0]):
            di,wi = d[i],w[i]
            wi[::5] = 0
            di *= wi
            dio,norm,res,m = O.calibration_omni.deconvolve_spectra2(di,wi,di.size/2)
            dio /= di.size
            _di = n.fft.ifft(di)
            _dio = n.fft.ifft(dio)
            _wi = n.fft.ifft(wi)
            _dic,info = a.deconv.clean(_di,_wi,tol=1e-6); _dic += info['res']
            #di.shape = (nchan,1)
            _di.shape = (nchan,1)
            wi.shape = (1,nchan)
            A = n.identity(nchan, dtype=n.complex)
            A = n.fft.ifft(n.fft.fft(A) * wi)
            p.subplot(111); C.arp.waterfall(A); p.show()
            At = A.T.conj()
            M = n.dot(At, A)
            M += alpha**2 * n.identity(M.shape[0], dtype=n.complex)
            _M = n.linalg.inv(M)
            _dim = n.dot(_M, n.dot(At,_di))
            p.subplot(121); C.arp.waterfall(M); p.colorbar()
            p.subplot(122); C.arp.waterfall(_M); p.colorbar()
            p.show()
            p.plot(n.real(_di)); p.plot(n.real(_dic)); p.plot(n.real(_dim)); p.plot(n.real(_dio)); p.show()
            
        _d = n.fft.ifft(d, axis=-1)
        p.subplot(121); C.arp.waterfall(d)
        p.subplot(122); C.arp.waterfall(_d)
        p.show()
