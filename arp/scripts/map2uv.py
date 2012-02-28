#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, ephem, os
import pylab as p
import time

o = optparse.OptionParser()
o.set_usage('map2img.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-o', '--output', dest='outfile', default='out.fits',
    help='Output filename.  Default out.fits')
o.add_option('--size', dest='size', type='int', default=300,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--mfreq', dest='mfreq', type='float', default=.150,
    help='Create an image evaluated at this frequency.  Default is .150 GHz')
opts, args = o.parse_args(sys.argv[1:])

JD = 2455746.2
ANT1 = 0
#ANT2s = [50, 60]
ANT2s = range(1, 10)
SDF0 = .1 / 1024
NCHAN = 1024 * 16
BL_SCALAR = 4.4

aa = a.cal.get_aa(opts.cal, .1/1024, .1, 1024)
afreqs = aa.get_afreqs()

print 'Reading', args[0],
t = time.time()
h = a.map.Map(fromfits=args[0])
h.set_interpol(True)
px = n.arange(h.npix())
tx,ty,tz = h.px2crd(px)
valid = n.where(tz > 0)
tx,ty,tz = tx[valid], ty[valid], tz[valid]
top = n.array([tx, ty, tz])
print time.time() - t
print '    NPIX:', h.npix()

if opts.pol != None:
    print '    Evaluating beam',
    t = time.time()
    aa.select_chans(n.array([512]))
    p1,p2 = opts.pol
    bm_resp = aa[ANT1].bm_response((tx,ty,tz), pol=p1) * n.conj(aa[ANT1].bm_response((tx,ty,tz), pol=p2))
    bm_resp = bm_resp.flatten()
    aa.select_chans()
    print time.time() - t
else: bm_resp = 1

phs_res, taubin = {}, {}
for ant2 in ANT2s:
    print '    Preparing', (ANT1,ant2),
    t = time.time()
    bl = a.miriad.ij2bl(ANT1,ant2)
    bx,by,bz = aa.get_baseline(ANT1, ant2, src='z')
    BL_LEN = n.sqrt(bx**2 + by**2 + bz**2) * BL_SCALAR
    SDF = 2**n.floor(n.log2(1 / BL_LEN / SDF0)) * SDF0
    BW = NCHAN * SDF
    bins = n.fft.fftfreq(NCHAN, SDF)
    bins = n.concatenate([bins[bins.size/2:], bins[:bins.size/2]])
    dbin = bins[1] - bins[0]
    tau = bx*tx + by*ty + bz*tz
    taubin[bl] = n.around(tau / dbin).astype(n.int)
    taures = tau - taubin[bl] * dbin
    phs_res[bl] = n.exp(-2j*n.pi*taures.astype(n.complex) * .150)
    taubin[bl].shape += (1,)
    print time.time() - t

for jd in n.arange(JD,JD+.01, .001):
    print jd
    aa.set_jultime(jd)

    print '    Precessing',
    t = time.time()
    m_precess = a.coord.convert_m('eq','eq', oepoch=aa.epoch)
    m = n.linalg.inv(n.dot(aa.eq2top_m, m_precess)) # takes top to J2000
    ex,ey,ez = n.dot(m, top)
    d = h[ex,ey,ez] * bm_resp
    print time.time() - t

    for ant2 in ANT2s:
        bl = a.miriad.ij2bl(ANT1,ant2)
        print ' ', (ANT1, ant2)
        print '    Misc preparations',
        t = time.time()
        bx,by,bz = aa.get_baseline(ANT1, ant2, src='z')
        BL_LEN = n.sqrt(bx**2 + by**2 + bz**2) * BL_SCALAR
        SDF = 2**n.floor(n.log2(1 / BL_LEN / SDF0)) * SDF0
        BW = NCHAN * SDF
        bins = n.fft.fftfreq(NCHAN, SDF)
        bins = n.concatenate([bins[bins.size/2:], bins[:bins.size/2]])
        dbin = bins[1] - bins[0]
        binedges = n.append(bins-dbin/2, bins[-1]+dbin/2)

        freq = n.fft.fftfreq(bins.size, bins[1] - bins[0])
        freq = n.where(freq < 0, freq + BW, freq)
        while n.all(freq < .15): freq += BW
        sdf = freq[1] - freq[0]
        window1 = n.where(freq + sdf/2 >= .05, 1., 0) * n.where(freq + sdf/2 < .25, 1., 0)
        #freq = freq.compress(window1)
        freq_pad = n.arange(.05,.25, .1/1024)
        sdf = freq_pad[1] - freq_pad[0]
        window2 = n.where(freq_pad + sdf/2 >= .1, 1., 0) * n.where(freq_pad + sdf/2 < .2, 1., 0)
        #freq_pad = freq_pad.compress(window2)
        print time.time() - t
        print '    Sub-delay bin phasing',
        t = time.time()
        d_bl = d * phs_res[bl]

        print time.time() - t
        print '    Binning',
        t = time.time()
        hist = n.zeros(bins.size, dtype=n.complex)
        a.utils.add2array(hist, taubin[bl], d_bl)
        #hist = n.histogram(tau, weights=d_bl, bins=binedges)[0]


        #p.subplot(122); p.semilogy(bins, n.abs(hist), 'gx')

        print time.time() - t
        print '    Computing coarse spectrum',
        t = time.time()
        hist = n.concatenate([hist[hist.size/2:], hist[:hist.size/2]])
        spec = n.fft.fft(hist).compress(window1)

        #p.subplot(121)
        #print spec.size
        #p.plot(freq, spec.real, 'g.')
        #p.xlim(.1,.2)

        #p.subplot(122)
        #dbins = n.fft.fftfreq(freq.size, freq[1]-freq[0])
        #p.semilogy(dbins, n.abs(dspec), 'g.')

        print time.time() - t
        print '    Computing fine spectrum',
        t = time.time()
        dspec = n.fft.ifft(spec)
        dspec_pad = n.zeros(2*afreqs.size, dtype=n.complex)
        dspec_pad[:dspec.size/2] = dspec[:dspec.size/2]
        dspec_pad[-dspec.size/2+1:] = dspec[-dspec.size/2+1:]
        spec_pad = n.fft.fft(dspec_pad)
        spec_pad = spec_pad.compress(window2)
        print time.time() - t

        #print spec_pad.shape

        #p.subplot(121)
        #p.plot(freq_pad, spec_pad.real, 'b-')

        #p.subplot(122)
        #dbins, dspec = n.fft.fftfreq(afreqs.size, afreqs[1]-afreqs[0]), n.fft.ifft(spec_pad)
        #p.semilogy(dbins, n.abs(dspec), 'b.')
#p.show()
