#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import capo as C
import sys

N = 5

t,d,f = C.arp.get_dict_of_uv_data(sys.argv[1:], '0_16', 'yy')

for bl in d:
  for pol in d[bl]:
    print d[bl][pol].shape
    U,s,V = n.linalg.svd(d[bl][pol], full_matrices=False)
    print s
    print U.shape, s.shape, V.shape

    p.subplot(N+1,3,1)
    p.subplot(N+1,3,2)
    C.arp.waterfall(d[bl][pol], drng=2)

    p.subplot(N+1,3,3)
    _d = C.arp.clean_transform(d[bl][pol], f[bl][pol], clean=1e-4)
    C.arp.waterfall(a.img.recenter(_d, (0, _d.shape[1]/2)), drng=2)
    
    for i in range(N):
        print i
        s_ = s.copy()
        s_[i+1:] *= 0
        S = n.diag(s_)
        dmdl = n.dot(U, n.dot(S, V))
        dres = n.where(f[bl][pol], 0, d[bl][pol] - dmdl)
    
        p.subplot(N+1,3,3*(i+1)+1)
        C.arp.waterfall(dmdl, drng=2)

        p.subplot(N+1,3,3*(i+1)+2)
        C.arp.waterfall(dres, drng=2)

        p.subplot(N+1,3,3*(i+1)+3)
        _d = C.arp.clean_transform(dres, f[bl][pol], clean=1e-4)
        C.arp.waterfall(a.img.recenter(_d, (0, _d.shape[1]/2)), drng=2)

    p.show()
#for i in u
