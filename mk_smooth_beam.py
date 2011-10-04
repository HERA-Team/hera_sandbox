#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, os, optparse
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_beam.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True,cal=True)
o.add_option('--tol', dest='tol', type='float', default=1.e-5,
    help="Tolerance for clean.  Default is 1e-5.")
o.add_option('--npix', dest='npix',type='int', default=1000,
    help="Number of pixels used in the sky image.")
o.add_option('--nmodes', dest='nmodes',type='int',default=6,
    help="Number of CLEAN components to use in the smooth reconstruction.")
o.add_option('--savefits', dest='savefits', action='store_true',
    help="Save the smooth beam as a HEALPix fits file.")
o.add_option('-m','--maxwgt',dest='maxwgt',type='float',default=1e6,
    help="Maximum weight in the deconvolution kernal.")
o.add_option('-o','--outfile',dest='outfile',type='string',default='output.npz',
    help="Name of the smooth beam coefficient npz file.")

opts,args = o.parse_args(sys.argv[1:])
center = opts.npix/2

print 'loading ', args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()

im = a.img.Img(size=(opts.npix*.4),res=.4)
x,y,z = im.get_top()

data = h.map[x,y,z]
kern = h.wgt[x,y,z]
oldkern = kern.copy()
kern = data.copy()
data *= (data/n.where(oldkern !=0,oldkern,1.))
old = data.shape
data.shape = x.shape
kern.shape = x.shape

maxwgt = n.max(kern)

mdata = n.where(x.mask == 1,0,data)
#mdata = n.where(kern < opts.maxwgt,mdata,mdata/n.where(kern!=0,kern,1)*opts.maxwgt)
#mkern = n.where(kern < opts.maxwgt,kern,opts.maxwgt)
mkern = n.where(x.mask == 1,maxwgt,kern)

mdata = a.img.recenter(mdata,(center,center))
mkern = a.img.recenter(mkern,(center,center))

_data = n.fft.fft2(mdata)
_kern = n.fft.fft2(mkern)

_beam = a.deconv.clean(_data,_kern,maxiter=10,stop_if_div=False,tol=opts.tol,verbose=True)

beam = n.fft.ifft2(_beam[0])*(_beam[0].size)
res = n.fft.ifft2(_beam[1]['res'])

_sbeam = a.img.recenter(_beam[0],(-center,-center))
_sbeam[:,:center-opts.nmodes] = 0.
_sbeam[:,center+opts.nmodes:] = 0.
_sbeam[:center-opts.nmodes,:] = 0.
_sbeam[center+opts.nmodes:,:] = 0.
#_sbeam[27:37,27:37] = 0.
sbeam = n.fft.ifft2(a.img.recenter(_sbeam,(center,center)))*(_sbeam.size)

outname = opts.outfile
#outname = args[0].replace('.fits','.smoothe.fits')
outfile = a.map.Map(h.nside())
xx,yy,zz = x.compressed(),y.compressed(),z.compressed()
#print xx.shape,yy.shape,zz.shape
hbeam = n.ma.array(a.img.recenter(sbeam,(-center,-center)),mask=x.mask).compressed()
hbeam.shape = xx.shape
outfile.add((xx,yy,zz),1,hbeam)
#norm = (outfile[0]+outfile[1]+outfile[2]+outfile[3])/4
#outfile.map.map/=norm

if opts.savefits == True:
    outfile.to_fits(outname,clobber=True)

coeffs = _sbeam[center-opts.nmodes:center+opts.nmodes,center-opts.nmodes:center+opts.nmodes]#/norm
#coeffs = _sbeam
print 'saving to ',outname.replace('.fits','.npz')
n.savez(outname.replace('.fits','.npz'),coeffs=coeffs)

ssbeam = n.fft.ifft2(a.img.recenter(_sbeam,(center,center)))*(_sbeam.size)


p.subplot(221)
p.imshow(n.abs(beam))
p.colorbar()
p.subplot(222)
p.imshow(n.real(res/(n.where(mkern!=0,mkern,1))))
p.colorbar()
p.subplot(223)
p.imshow(n.abs(ssbeam),interpolation='Nearest')
p.colorbar()
p.subplot(224)
p.imshow(mkern)
p.colorbar()

p.show()
