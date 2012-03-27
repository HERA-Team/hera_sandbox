#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, sys

times, d, f = C.arp.get_dict_of_uv_data(sys.argv[1:], '0_16', 'xx', verbose=True)

for bl in d:
    for cnt in range(len(times)):
        print cnt
        C.arp.sinuspike(d[bl][cnt], f=f[bl][cnt])
        if cnt > 5: break
