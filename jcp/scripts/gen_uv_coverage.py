#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts, args = o.parse_args(sys.argv[1:])

def beamgridder(sigma,xcen,ycen,size):
    #makes FFT of 2D gaussian beam with sigma; returns an array of size x size
    crds = n.mgrid[0:size,0:size]
    cen = size/2 - 0.5 # correction for centering
    xcen += cen
    ycen = -1*ycen + cen
    beam = n.exp(-((crds[1]-xcen)/(2*sigma))**2)*n.exp(-((crds[0]-ycen)/(2*sigma))**2)
    return beam/n.sum(beam)

#observing time
cen_jd = 2454600.90911
obs_duration = 40. #minutes
start_jd = cen_jd - (1./24)*((obs_duration/60)/2)
end_jd = cen_jd + (1./24)*(((obs_duration-1)/60)/2)
times = n.arange(start_jd,end_jd,(1./24/60))

#other observing parameters
fq = .150
dish_size_in_lambda = 5.51 #wavelengths; PAPER primary beam model for HERA
SIZE = 600 #uv plane size in wavelengths
#sigma = 0.0727
#sigma = dish_size_in_lambda/2.35 #in uv plane
sigma = 1/(2*n.pi*0.0727) #analytic with right Fourier convention

#load array information
aa = a.cal.get_aa(opts.cal,n.array([.150]))
cat = a.src.get_catalog(opts.cal,'z')
nants = len(aa)

dim = n.round(SIZE/dish_size_in_lambda/2)*2 - 1 # round to nearest odd
uv = n.zeros((dim,dim))
print 'uv plane is (%i,%i) pixels' % (dim,dim)

aa.set_jultime(cen_jd)
obs_lst = aa.sidereal_time()
obs_zen = a.phs.RadioFixedBody(obs_lst,aa.lat)
obs_zen.compute(aa)

mags = []
for cnt, t in enumerate(times):
    print 'Working on integration %i of %i' % (cnt+1, obs_duration)
    for i in xrange(nants):
        #print i
        for j in xrange(nants):
            if i == j: continue #don't bother with autocorrelations
            #if i > j: continue #don't double count
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            obs_zen.compute(aa)
            u,v,w = aa.gen_uvw(i,j,src=obs_zen)
            _beam = beamgridder(sigma=sigma,xcen=u/dish_size_in_lambda,ycen=v/dish_size_in_lambda,size=dim)
            #print sigma, u/dish_size_in_lambda,v/dish_size_in_lambda, dim
            uv += _beam

#uv[:,:dim/2] = 0
#uv[dim/2:,dim/2] = 0

n.save('uvcov.npy',uv)

print 'there are %i minutes of integration in the uv plane' % n.sum(uv)

p.imshow(uv,interpolation='nearest')
p.colorbar()
p.show()

