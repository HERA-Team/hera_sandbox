#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, scipy
from mpl_toolkits.basemap import Basemap

FQ = .151
#BINWIDTH = 1. / 7200. # Hz
BINWIDTH = 1. / 3600. # Hz
bin_edges = n.arange(-.01+BINWIDTH/2,.01,BINWIDTH)
NBINS = len(bin_edges)
DIM = 400
RES = .5
SIZE = DIM / RES

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))

h = a.healpix.HealpixMap(nside=64)
im = a.img.Img(DIM, res=RES)

def proj_hmap(dat, mode='e'):
    h.map = dat
    if mode == 'e': x,y,z = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, 0)
    return d

m = Basemap(projection='ortho', lat_0=aa.lat, lon_0=180, rsphere=1.)

def plot_hmap(dat, mode='log', mx=0, drng=3):
    d = proj_hmap(dat,'e')
    if mode == 'log': d = n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
    elif mode == 'real': d = d.real
    m.imshow(d, vmax=mx, vmin=mx-drng, origin='lower', interpolation='bicubic', cmap='jet')

def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / h.npix()

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmI = 0.5 * (_bmx**2 + _bmy**2)
bmI = n.where(tz > 0, bmI, 0)

#def nos(dat, t):
#    pwr_bm = beam_area(bm)
#    pwr2_bm = beam_area(dat**2)
#    return pwr_bm**2 / pwr2_bm / t

#print beam_area(bm), beam_area(bm**2), nos(bm, 1.), nos(bm,1.) / n.sqrt(NBINS)

xyz = (xyz[1],xyz[0],xyz[2])
#bl = n.array([100, 0, 0])
bl = aa.get_baseline(0,16,'r') * FQ
#bl = aa.get_baseline(0,28,'r') * FQ
#bl = aa.get_baseline(0,3,'r') * 10 * FQ
print 'Baseline:', bl
fng = C.frf_conv.mk_fng(bl, *xyz)
frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,16))
bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,16), ref_chan=0)
tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)
skypass = n.where(n.abs(bins)<0.5/42.9, 1, 0)
crosstalk = n.ones_like(frbins)
#crosstalk[14:] = 0 # 600s / 42.9s = 14
#crosstalk[3*14:] = 0 # 600s / 42.9s = 14
#crosstalk[:3*14] *= a.dsp.gen_window(3*14, 'blackman-harris')
#crosstalk[3*84:] = 0 # 3600s / 42.9s = 84
#crosstalk[:3*84] *= a.dsp.gen_window(3*84, 'blackman-harris')
crosstalk = n.fft.fftshift(n.abs(n.fft.ifft(crosstalk)))
crosstalk /= crosstalk.max()
print frfs[0].max()
#frfs[0] /= frfs[0].max() # XXX shouldn't have to do this
p.plot(frbins,1-crosstalk)
p.plot(bins,frf)
p.plot(frbins,frfs[0])
#frfs[0] *= (1-crosstalk)
p.show()
wgt = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear')
fng_wgt = wgt(fng)
fng_bm = fng_wgt * bmI
#print fng.max(), fng.min()
noise_ratio = n.sum(n.abs(frfs[0])**2)/n.sum(n.abs(skypass)**2)
print noise_ratio, 42.9 / noise_ratio, n.sqrt(noise_ratio)


#hist, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmI**2)
#hist, bin_edges = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm**2) # or bm**2?
f = 0
for cnt,b1 in enumerate(bin_edges[:-1]):
    b2 = bin_edges[cnt+1]
    if cnt % 2 == 0: f = n.where(n.logical_and(fng >= b1, fng < b2), 1., f)

p.subplot(131)
plot_hmap(bmI)

p.subplot(132)
#p.imshow(proj_hmap(fng,'e'), vmax=0, vmin=-.007, origin='lower')
plot_hmap(fng_wgt)
#
p.subplot(133)
plot_hmap(fng_bm)
#
#p.subplot(234)
##hist, bins = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm) # or bm**2?
#bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#hist1 = hist/n.sum(hist)
#print 'Hist:', hist1.max()
##p.plot(bins, hist1/hist1.max())
#p.plot(bins, hist1/hist1.max())
#p.plot(bins, n.sqrt(hist1/hist1.max()))
#
#noise_sum, noise_wgt = 0, 0
#wgts = []
#f = 0
#for cnt,b1 in enumerate(bin_edges[:-1]):
#    b2 = bin_edges[cnt+1]
#    fng_wgt = n.where(n.logical_and(fng >= b1, fng < b2), 1., 0)
#    fng_amp = bm * fng_wgt
#    _n = nos(fng_amp,NBINS)
#    w = 1. / _n**2
#    wgts.append(w)
#    if w > 0:
#        #print cnt, b1, _n, w
#        noise_sum += (_n*w)**2
#        noise_wgt += w
#    f = n.where(n.logical_and(fng >= b1, fng < b2), w, f)
#wgts = n.array(wgts)
#
#p.plot(bins, wgts/wgts.max())
#
#fng_amp = bm * f
#print fng_amp.max()
#fng_amp /= fng_amp.max()
bmI_area = beam_area(bmI**2)
fng_bm_area = beam_area(fng_bm**2)
print 'Whole PP beam area:', bmI_area
print 'Filtered PP beam area:', fng_bm_area
print 'Correction factor:', bmI_area / fng_bm_area
#print n.sqrt(noise_sum) / noise_wgt
#
#print 'Improvement:', (n.sqrt(noise_sum)/noise_wgt) / (nos(bm,1.) / n.sqrt(NBINS))
#    
#p.subplot(235)
#p.imshow(proj_hmap(f,'e'), origin='lower')
#
#p.subplot(236)
#p.imshow(proj_hmap(fng_amp,'e'), origin='lower')

p.show()
