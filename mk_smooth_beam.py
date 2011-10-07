#!/usr/bin/env python

import aipy as a, numpy as n, sys, optparse

o = optparse.OptionParser()
o.add_option('--tol', dest='tol', type='float', default=1.e-5,
    help="Tolerance for clean.  Default is 1e-5.")
o.add_option('--nmodes', dest='nmodes',type='int',default=6,
    help="Number of CLEAN components to use in the smooth reconstruction.")
o.add_option('--savefits', dest='savefits', action='store_true',
    help="Save the smooth beam as a HEALPix fits file.")
o.add_option('-o','--outfile',dest='outfile',type='string',default='output.npz',
    help="Name of the smooth beam coefficient npz file.")
opts,args = o.parse_args(sys.argv[1:])
filename = args[0]

print 'Reading:', filename
h = a.map.Map(fromfits=filename)
im = a.img.Img(size=400,res=.4)
x,y,z = im.get_top()

data,kern = h.map[x,y,z], h.wgt[x,y,z]
data.shape, kern.shape = x.shape, x.shape
data = n.where(x.mask, 0, data)
# Reweight data for optimal SNR
kern,data = data.copy(), data**2 / n.where(kern == 0, 1, kern)
kern = n.where(x.mask, kern.max(), kern)

_data,_kern = n.fft.fft2(data), n.fft.fft2(kern)
# Only allow clean components out to specified number of modes
#area = n.ones(_data.shape, dtype=n.int)
#area[opts.nmodes+1:-opts.nmodes] = 0
#area[:,opts.nmodes+1:-opts.nmodes] = 0
_beam,info = a.deconv.clean(_data, _kern, #area=area,
    maxiter=10, stop_if_div=False, tol=opts.tol, verbose=True, gain=.8)
# Only allow clean components out to specified number of modes
_beam[opts.nmodes+1:-opts.nmodes] = 0
_beam[:,opts.nmodes+1:-opts.nmodes] = 0
beam = n.fft.ifft2(_beam)*_beam.size
beam = n.where(x.mask, 0, beam)
res = data - kern * beam

if False:
    import pylab as p
    for i,d in enumerate([data/n.where(kern == 0, 1, kern),kern,beam,res/n.where(kern == 0, 1, kern)]):
        p.subplot(2,2,i+1)
        d = a.img.recenter(d,(500,500))
        d = n.log10(n.abs(d))
        p.imshow(d, vmax=d.max(), vmin=d.max()-2)
    p.show()

outname = filename.replace('.fits','_sm.fits')
npzname = filename.replace('.fits','_sm.npz')
outfile = a.map.Map(h.nside())
hbeam = n.ma.array(beam,mask=x.mask).compressed()
x,y,z = x.compressed(),y.compressed(),z.compressed()
outfile.add((x,y,z),1,hbeam)

print 'Saving smooth beam to', outname
outfile.to_fits(outname,clobber=True)

center = n.array(_beam.shape)/2
_beam = a.img.recenter(_beam,center)
coeffs = _beam[center[0]-opts.nmodes:center[0]+opts.nmodes+1,center[1]-opts.nmodes:center[1]+opts.nmodes+1]

print 'Saving coefficients to', npzname
n.savez(npzname, coeffs=coeffs)

