#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p

def mk_fng(bl, ex, ey, ez):
    #f = -(bl[0]*n.sin(lons) + bl[1]*n.cos(lons)) * n.cos(lats)
    return -2*n.pi/a.const.sidereal_day * (bl[0]*ex + bl[1]*ey) * n.sqrt(1 - ez**2)

aa = a.cal.get_aa('psa898_v002', n.array([.150]))

DIM = 400
RES = .5
SIZE = DIM / RES
im = a.img.Img(DIM, res=RES)
tx,ty,tz = im.get_top(center=(SIZE/2,SIZE/2))
SH = tx.shape
valid = n.logical_not(tx.mask)
tx,ty,tz = tx.flatten(), ty.flatten(), tz.flatten()
amp = aa[0].bm_response((tx,ty,tz),pol='x')**2
amp.shape = SH
amp = n.where(valid, amp, 0)
amp /= amp.sum()

p.subplot(231)
p.imshow(amp, origin='lower')

ex,ey,ez = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
#ex,ey,ez = ex.flatten(), ey.flatten(), ez.flatten()

bl = n.array([100, 0, 0])

p.subplot(232)
fng = mk_fng(bl, ex, ey, ez)
fng = n.where(valid, fng, 0)
print fng.max(), fng.min()
p.imshow(fng, vmax=0, vmin=-.007, origin='lower')

p.subplot(233)
f = (fng / .001) % 2
f = n.where(f > 1, 1, 0)
p.imshow(f*amp, origin='lower')

p.subplot(234)
NBINS = 200
hist, bins = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=amp)
bins = 0.5 * (bins[:-1] + bins[1:])
hist1 = hist/n.sum(hist)
p.plot(bins, hist1/hist1.max())
hist2 = hist**2
hist2 /= hist2.max()
p.plot(bins, hist2)
all_sky = n.where(hist1 > 0, 1., 0)
#hist2 = n.where(hist2 > 0.5 * hist2.max(), 1., 0)
frac_int_time = all_sky.sum() / hist2.sum()
#p.plot(bins, n.cumsum(hist2))

import scipy.interpolate
wgts = scipy.interpolate.interp1d(bins, hist2/hist2.max(), kind='cubic')

fng_wgt = wgts(fng)

p.subplot(235)
p.imshow(fng_wgt, origin='lower')

p.subplot(236)
fng_amp = amp * fng_wgt
#fng_amp /= fng_amp.sum()

if False:
    bm_area = n.sqrt(n.sum(amp**2))
    fng_bm_area = n.sqrt(n.sum(fng_amp**2))
else:
    bm_area = n.sum(amp)
    fng_bm_area = n.sum(fng_amp)
frac_fng_bm_area = fng_bm_area / bm_area

print 'Reduction in beam area:', frac_fng_bm_area
print 'Increase in integration time:', frac_int_time
print 'Net SNR improvement:', frac_fng_bm_area * frac_int_time

p.imshow(fng_amp, origin='lower')

p.show()
