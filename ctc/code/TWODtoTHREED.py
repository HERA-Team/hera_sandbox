## NAME: 
#         TWODtoTHREED.py
## PURPOSE: 
#         Takes a 2D image and creates a 3D healpix map

#------------------------------------------------------------------------------------------------

import aipy
import numpy
import pylab
import pyfits

"""

#make 2D image

size = 100.0
pic = numpy.random.random((size,size)) #size is 100x100
pylab.imshow(pic)
#pylab.show()

"""

#take 2D image

image = pyfits.open("/Users/carinacheng/Desktop/Carina/UCBResearch/images/modelpsf_ghezgroup.fits")
pic = image[0].data
size = len(pic[0])

#give 3D coordinates to every pixel
x = numpy.arange(-1.0,1.0, 1.0/(size)) #1D array of length size that goes from -1 to 1
x = numpy.resize(x,(size,size)) #resizes array to be 2D
y = numpy.arange(-1.0,1.0, 1.0/(size))
y = numpy.resize(y,(size,size))
y = y.transpose() #transposes array so that it goes through the columns instead of rows
z = 1-x**2-y**2
mask = numpy.where(z<0, 1, 0) #masks out negative z's (puts 1's there and 0's otherwise)
z = numpy.ma.array(numpy.sqrt(z), mask=mask) #masked array function where any 'true' value (1's) are excluded in computations
x = numpy.ma.array(x, mask=mask)
y = numpy.ma.array(y, mask=mask)
x = x.compressed() #gets rid of masked values because the C code later on doesn't see them
y = y.compressed()
z = z.compressed() 

#make healpix 3D map

bm = aipy.map.Map(nside=512) #bm is now a Map class object, nside=# of pixel regions on sphere (must be power of 2)
   #note: bm.map is a healpix map (type bm.map? to confirm)
#px = bm.crd2px(x,y,z) #not necessary here actually (equivalent to bm.add command)
flx = numpy.ma.array(pic, mask=mask).compressed() #1D array of the original 2D image, ignoring masked values
wgt = numpy.ones_like(flx) #1D array of 1's for the weight of the same size/type as flx array
bm.add((x,y,z), wgt, flx) #can use 'add' or 'put' (pretty much equivalent)

bm.to_fits('/Users/carinacheng/Desktop/Carina/UCBResearch/images/3Dfg.fits', clobber=True) #clobber overwrites the fits file each time

   #type 'plot_map.py 3Dfg' in terminal to plot it (--nobar if colorbar error)
