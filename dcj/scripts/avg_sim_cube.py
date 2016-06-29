#! /usr/bin/env  python

from dcj_etc import downsample
from cosmo_units import *
import optparse,sys,time
from scipy import ndimage

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--smooth', default=10.,type=float,
    help='Smooth in the image plane by this many arcminutes wide gaussian. [default 10]')
opts, args = o.parse_args(sys.argv[1:])
flush = sys.stdout.flush
Dx = 1200. #proper distance size of a cube
N = 1024
z_end = 7.1 #ending redshift of the cube
print "loading cosmology with O_l = %f, O_M = %f, Ho = %f"%(O_l,O_M,Ho)
zi = n.linspace(6,15)
D_i = [DM(z) for z in zi]
dtheta = n.poly1d(n.polyfit(D_i,[n.degrees(r2theta(Dx/1024,z))*60 for z in zi],4)) #pixel size in arcminutes
z_p = n.polyfit(D_i,zi,4)
z = n.poly1d(z_p)

dt = n.degrees(r2theta(Dx/1024,10))*60
print dt
pix_psf = n.ceil(opts.smooth/dt)
pix_resample = int(max([pix_psf/4,1]))
print "smoothing to %4.2f arcminutes. (%d pixels)"%(opts.smooth,pix_psf)
print "resampling image plane to %d pixels(%5.3f arcmin) "%(pix_resample, pix_resample*dt)
flush()
print "making empty downsampled array",;flush()
T_obs = n.zeros((N/pix_resample,N/pix_resample,N))
print "done"

print "reading T file...",
flush()
T = n.fromfile('1024_1200Mpc_lc_z7.1.temp',dtype='<f8')
T.shape = (N,N,N)
print "done"

if False:
    print "computing cosmology",
    R = DM(z_end) + n.arange(N)*Dx/N  #the proper distances of the slices
    Z = z(R) #the redshifts of the slices
    f = f21/(1+Z) #the frequencies of the slices
    print "done"
#print "computing pixel size vs redshift"; flush()
#t0 = time.time()
#print dtheta(R[0])
#print "this will take %5.3f minutes "%((time.time() - t0)/60*N);flush()
#dthetas = n.degrees([r2theta(R[i],Z[i]) for i in range(N)])*60 #pixel size in arcminutes of the slices
#dthetas = n.array([dtheta(r) for r in R])
#print "done";flush()
print "smoothing data slices"
for i in range(N):
    temp = ndimage.gaussian_filter(T[:,:,i],pix_psf)
    T_obs[:,:,i] = downsample(temp,pix_resample)
    if i%64==0: print i/float(N),
    else: print '.',
    flush()
T_obs.tofile('1024_1200Mpc_lc_z7.1.temp_smoothsmall')

