#! /usr/bin/env python
import aipy as a, capo as C, numpy as n, pylab as p
import sys, pyfits

class HealpixMap(a.healpix.HealpixMap):
    '''Subclass to match the slightly different Healpix format James wrote.'''
    def from_fits(self, filename, hdunum=1, colnum=0):
        hdu = pyfits.open(filename)[hdunum]
        data = hdu.data.field(colnum)
        data = data.flatten()
        if not data.dtype.isnative:
            data.dtype = data.dtype.newbyteorder()
            data.byteswap(True)
        scheme= hdu.header['ORDERING'][:4]
        self.set_map(data, scheme=scheme)

BL = 14e2 # cm
bl_ns = BL / a.const.len_ns
bl = BL * n.array([1.,0,0])

h = HealpixMap(fromfits=sys.argv[-1])
x,y,z = h.px2crd(n.arange(h.npix()))
bm = h[x,y,z]
#bl_mag = n.sqrt(BL**2 - (bl[0]*x+bl[1]*y+bl[2]*z)**2)
bl_mag = n.sqrt(BL**2 - (bl[0]*z+bl[1]*y+bl[2]*x)**2)


h.map = bl_mag; C.arp.plot_hmap_ortho(h, mode='real'); p.show()
Tsync = 230. * ((bl_mag+3) / 2e2)**-3 # K
Tpnt = 10 # K
bm = n.where(z > 0, bm, 0)
hist,bins = n.histogram(x, weights=(Tsync+Tpnt)*bm, bins=100)
hist /= hist.max()
bins = .5*(bins[:-1]+bins[1:])
taus = bins * bl_ns


p.plot(taus, hist)
p.show()
