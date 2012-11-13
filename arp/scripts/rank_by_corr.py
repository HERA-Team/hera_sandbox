#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p
import sys

def corr(d1, d2):
    c = n.abs(n.average(d1*n.conj(d2)))**2
    c /= n.average(n.abs(d1) * n.abs(d2))**2
    return n.sqrt(c)
       

#p.ion()
SZ = 64*65/2

grps = {}
grp_dat = {}

data = {}
for filename in sys.argv[1:]:
    print filename
    uv = a.miriad.UV(filename)
    bls,curtime = {},None
    for (crd,t,(i,j)),d1,f1 in uv.all(raw=True):
        if a.miriad.pol2str[uv['pol']] != 'xx': continue
        if t != curtime:
            if not curtime is None: break
            bls,curtime = {},t
            print t
        bl1 = a.miriad.ij2bl(i,j)
        data[bl1] = d1
        #px1 = len(bls)
        #for px2,(bl2,d2) in enumerate(bls.items()):
        mx_grp,mx_scr = None, 0.
        for g in grps:
            bl2 = grps[g][0]
            d2 = grp_dat[g]
            c = corr(d1, d2)
            if c > mx_scr: mx_scr,mx_grp = c,g
        if mx_scr > .5:
            grps[mx_grp].append(bl1)
            print '%2d-%2d/%2d: %5.2f' % (a.miriad.bl2ij(bl1) + (mx_grp,mx_scr))
        else:
            g = len(grps)
            grps[g] = [bl1]
            grp_dat[g] = d1
    C = n.zeros((SZ,SZ), dtype=n.float)
    offset = 0
    for g in grps:
        print g
        bls = []
        for bl1 in grps[g]:
            for px2,bl2 in enumerate(bls):
                px1 = offset+len(bls)
                px2 += offset
                C[px1,px2] = corr(data[bl1], data[bl2])
                print '%2d-%2d/%2d-%2d: %5.2f' % (a.miriad.bl2ij(bl1) + a.miriad.bl2ij(bl2) + (C[px1,px2],))
            bls.append(bl1)
        offset += len(bls)
            
    #plt = p.imshow(C, vmax=1, vmin=0)
    p.imshow(C, vmax=1, vmin=0)
    p.show()

    
        #if px1 % 100 == 0:
        #    print px1
        #    plt.set_data(C)
        #    p.draw()
        #bls[bl1] = d1
