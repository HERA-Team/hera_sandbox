#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, os, optparse
from mpl_toolkits.basemap import Basemap
import beamuv

o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True)
o.add_option('--nmodes', dest='nmodes',type='int',default=6,
    help="Number of CLEAN components to use in the smooth reconstruction.")
o.add_option('-s','--size',dest='size',type='int',default=1000,
    help="Size of the uv matrix in pixels.  Default is 1000.")

opts,args = o.parse_args(sys.argv[1:])
center = opts.size/2

_coeffs = beamuv.coeffs_from_file(args[0])
offset = _coeffs.shape[0]/2
_beam = n.zeros((opts.size,opts.size),dtype=n.complex64)
_beam[center-offset:center+offset+1,center-offset:center+offset+1] = _coeffs

_sbeam = _beam.copy()
_sbeam[:,:center-opts.nmodes] = 0.
_sbeam[:,center+opts.nmodes+1:] = 0.
_sbeam[:center-opts.nmodes,:] = 0.
_sbeam[center+opts.nmodes+1:,:] = 0.

beam = n.fft.ifft2(a.img.recenter(_beam,(-center,-center)))*_beam.size
sbeam = n.fft.ifft2(a.img.recenter(_sbeam,(-center,-center)))*_sbeam.size

outnamea,outnameb = args[0].split('_')
outname = outnamea + '_f%s_' % opts.nmodes + outnameb
print 'saving smooth beam coeffecients to', outname

coeffs = _sbeam[center-opts.nmodes:center+opts.nmodes+1,center-opts.nmodes:center+opts.nmodes+1]
n.savez(outname,coeffs=coeffs)

if True:
    outname = outname.replace('.npz','.fits')
    outfile = a.map.Map(nside=64)
    im = a.img.Img(size=400,res=.4)
    x,y,z = im.get_top()
    hbeam = n.ma.array(sbeam,mask=x.mask).compressed()
    #hbeam = n.ma.array(a.img.recenter(sbeam,(-center,-center)),mask=x.mask).compressed()
    x,y,z = x.compressed(),y.compressed(),z.compressed()
    outfile.add((x,y,z),1,hbeam)
    outfile.to_fits(outname,clobber=True)

p.subplot(121)
p.imshow(n.abs(beam))
p.colorbar()
p.subplot(122)
p.imshow(n.abs(sbeam))
p.colorbar()

p.show()
