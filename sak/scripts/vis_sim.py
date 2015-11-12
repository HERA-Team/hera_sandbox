import pylab, random, os, sys
import healpy as hp
import numpy as np
from astropy import units as u
from astropy import constants as c
from scipy import integrate
from bm_prms import prms

assert(sys.argv[1] != None)

def rotate_hmap(map,rot):
	npix = map.shape[0]
	nside = hp.npix2nside(npix)

	rotmap = np.zeros(npix)
	ipix = np.arange(npix)
	t,p = hp.pix2ang(nside,ipix)

	r = hp.Rotator(rot=rot)

	# For each pixel in the new map, find where it would have come 
	# from in the old    
	trot,prot = r(t,p)
	ipix_rot = hp.ang2pix(nside,trot,prot)

	rotmap = map[ipix_rot]

	return rotmap


freqs = np.linspace(0.117,0.182,num=131) #aipy likes GHz units. avoiding band edges
nside = 128
npix = hp.nside2npix(nside)

pwd = os.getcwd()
if os.path.exists(pwd+'/XX_beam_maps.npz'):
	beam = np.load('XX_beam_maps.npz')['maps']
	
else:
	#create beams 
	beams = np.zeros((freqs.shape[0],npix))

	###CLEARLY THIS ONLY WORKS FOR LINEAR POLARIZATIONS XX AND YY

	print 'Calculating beams:'
	for i, freq in enumerate(freqs):
		print freq,'GHz'
		bm = prms['beam'](np.array([freq]),nside=nside,lmax=20,mmax=20,deg=7)
		bm.set_params(prms['bm_prms'])
		px = range(hp.nside2npix(nside))
		xyz = hp.pix2vec(nside,px)
		poly = np.array([h.map[px] for h in bm.hmap])
		Axx = np.polyval(poly,freq)
		Axx = np.where(xyz[-1] >= 0, Axx, 0)
		Axx /= Axx.max()
		Axx = Axx*Axx
		beams[i,:] = rotate_hmap(Axx,[0,90]) #[0,0]=north pole, [0,90]=equator

	np.savez('XX_beam_maps.npz',maps=beams)

	beam = beams



#calculate relevant map parameters
c = 3e8 #m/s
#nside = hp.npix2nside(beam.shape[0])
#npix = beam.shape[0]
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

#we care about degree scales ~21 degrees
#lmax=9
lmax=3*nside - 1
l,m = hp.Alm.getlm(lmax)

#frequencies
nfreq=freqs.shape[0]
nu = np.outer(np.linspace(100e6,200e6,num=nfreq),np.ones(npix))#*u.Hz


#define sky -- completely arbitrary choice of temp
uniform_sky = np.ones(npix)*100.#*u.K

#completely arbitrary choice of noise level
noise = np.zeros(npix)
for i in range(npix): noise[i] = random.uniform(-100,100)
#noise = noise * u.K

####uniform sky tests and point source tests

#sky = uniform_sky+noise
#sky = uniform_sky
"""
###DEBUG
###sky = np.zeros(npix)

#define a point source ###in both hemispheres
theta0 = np.pi/2.
phi0 = 0.
pix0 = hp.ang2pix(nside,theta0,phi0)
nbs = hp.get_all_neighbours(nside,theta0,phi=phi0)
sky[pix0]=500#*u.K
for nb in nbs: sky[nb]=500#*u.K
sky = rotate_hmap(sky,[180,0])


#for nb in nbs: sky[nb]=500#*u.K

theta0 = np.pi/2.
phi0 = 0.
pix0 = hp.ang2pix(nside,theta0,phi0)
nbs = hp.get_all_neighbours(nside,theta0,phi=phi0)
sky[pix0]=1000#*u.K
"""

########Verifying Ridhima sims
sky = hp.read_map('stokesI-f100_j2455819.54472.fits')


#hp.orthview(sky)
#pylab.show()

#promote sky to matrix for frequency axis
sky = np.outer(np.ones(nfreq),sky)*pow(nu/150e6,-0.7)

#decompose sky into alm's
n_alm = len(m)
alm = np.zeros((nfreq,n_alm),dtype='complex128')
print 'Calculating sky a_lm values:'
for i in range(nfreq):
	print nu[i,0]/1e6,'MHz'
	alm[i,:] = hp.map2alm(sky[i,:],lmax=lmax,iter=3)

#calculate fringe factor (true for all freqs)	
s  = np.array(hp.pix2vec(nside,ipix))

#
#
bl_length = 100.
#
#

#b = np.resize(np.repeat(np.array([0,bl_length,0]),npix),[3,npix])#*u.meter

#Checking Ridhima's simulations
b = np.resize(np.repeat(np.array([27.322,-239.994,0]),npix),[3,npix])#*u.meter


b_dot_s = np.sum(b*s,axis=0)
factor = np.exp(1.j*np.outer(np.ones(nfreq),b_dot_s)*nu/c) #c.c didn't work

#decompose B*fringe into alm's
blm = np.zeros((nfreq,n_alm),dtype='complex128')
print 'Calculating instrument a_lm values:'
for i in range(nfreq):
	print nu[i,0]/1e6,'MHz'
	blm[i,:] = hp.map2alm((beam*factor)[i,:],lmax=lmax,iter=3)

#phasing
#rot_ang = np.linspace(0,5.0592,num=1632)<--LST bin rate
rot_ang = np.linspace(-np.pi,np.pi,num=360*4) #<--24 hours
n = len(rot_ang)

vis = np.zeros([n,nfreq],dtype='complex128')

print 'Calculating visibilities. Time stamp:'
for i in range(n):
	print i
	rotation = np.outer(np.ones(nfreq),np.exp(-1.j*m*rot_ang[i]))
	vis[i,:] = np.sum(alm*blm*rotation,axis=1)

savefile = sys.argv[1]
print 'Saving visibility to %s...'%savefile
np.savez(savefile,vis=vis)	
