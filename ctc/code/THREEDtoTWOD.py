## NAME: 
#         THREEDtoTWOD.py
## PURPOSE: 
#         Takes a 3D healpix map and plots the 2D projection centered on a particular RA and dec
#         Interactive mode: 2D image rotates (RA changes)

#------------------------------------------------------------------------------------------------

import aipy
import numpy
import pylab
import pyfits

pylab.ion() #interactive mode on

plt = None #just a variable

for ii in numpy.arange(0,2*numpy.pi,0.1): #for loop for ra from 0 to 2pi with steps of 0.1

    #get 3D coordinates for every pixel in a 2D image
    
    size = 300 #number of wavelengths when starting in uv domain
    
    img2d = aipy.img.Img(size=size, res = .5)
    crd = img2d.get_eq(ra=ii,dec=-(numpy.pi)/12.0) #equatorial coordinates (bug at ra=pi!)
    x,y,z = crd
    
    sh = x.shape #remember dimensions of x
    mask = x.mask #remember where the masked values are
    
    #x,y,z = x.filled(0), y.filled(0), z.filled(1) #need to put good coordinates in, but these get masked out later anyways with fluxes
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    
    #get array of values for image
    
    img3d = aipy.map.Map(fromfits = '/Users/carinacheng/Desktop/Carina/UCBResearch/images/lambda_haslam408_dsds_eq.fits') #reads in 3D image
    
    fluxes = img3d[x,y,z] #gets fluxes already weighted (1D array)
    #fluxes = img3d.get(crd)[1] #gets fluxes at the coordinates but they aren't weighted
    
    #plot 2D image
    
    fluxes.shape = sh
    fluxes = numpy.ma.array(fluxes, mask=mask)

    if plt == None: #for the first time, plot as specified below

        plt = pylab.imshow(numpy.fft.fftshift(numpy.log10(fluxes)),interpolation='nearest',origin='lower',extent=(-1,1,-1,1)) 
        pylab.show()

    else: #update the plot (redraw it) at each step
        
        plt.set_data(numpy.fft.fftshift(numpy.log10(fluxes)))
        pylab.draw()






