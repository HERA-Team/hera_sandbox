#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts, args = o.parse_args(sys.argv[1:])

array = 'baobab'

if array == 'baobab':
    obs_duration = 30. # minutes
    dish_size_in_lambda = 5 #simple lambda/D
    SIZE = 300 #uv plane size in wavelengths
    sigma = 1. #not used with delta function gridding

if array == 'fidrsd':
    obs_duration = 80. # minutes
    dish_size_in_lambda = 3 #simple lambda/D
    SIZE = 600 #uv plane size in wavelengths
    sigma = 1. #not used with delta function gridding

if array == 'hera':
    obs_duration = 40. # minutes
    dish_size_in_lambda = 7 #simple lambda/D
    SIZE = 600 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.0643)

if array == 'paper':
    obs_duration = 120. # minutes
    dish_size_in_lambda = 1.5 #simple lambda/D #you published a paper with 2
    SIZE = 300 #uv plane size in wavelengths
    #SIZE = 50 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.2251)

if array == 'mwa':
    obs_duration = 120. # minutes
    dish_size_in_lambda = 2.65 #simple lambda/D
    SIZE = 600 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.1702)

if array == 'lofar':
    obs_duration = 15. # minutes
    dish_size_in_lambda = 15.4 #simple lambda/D
    SIZE = 600 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.0293)

if array == 'skalo1':
    obs_duration = 15. # minutes
    dish_size_in_lambda = 17.5 #simple lambda/D
    SIZE = 600 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.02577)

if array == 'skamid1':
    obs_duration = 40. # minutes
    dish_size_in_lambda = 7.5 #simple lambda/D
    SIZE = 2000 #uv plane size in wavelengths
    sigma = 1/(2*n.pi*0.02577)

def beamgridder(sigma,xcen,ycen,size):
    #makes FFT of 2D gaussian beam with sigma; returns an array of size x size
    crds = n.mgrid[0:size,0:size]
    cen = size/2 - 0.5 # correction for centering
    xcen += cen
    ycen = -1*ycen + cen
    #beam = n.exp(-((crds[1]-xcen)**2/(2*sigma**2)))*n.exp(-((crds[0]-ycen)**2/(2*sigma**2))) #beam gridder
    beam = n.zeros((size,size))
    if n.abs(round(ycen)) > size - 1 or n.abs(round(xcen)) > size - 1: 
        return beam
    else:
        beam[round(ycen),round(xcen)] = 1. #single pixel gridder
        return beam/n.sum(beam)

#observing time
cen_jd = 2454600.90911
start_jd = cen_jd - (1./24)*((obs_duration/60)/2)
end_jd = cen_jd + (1./24)*(((obs_duration-1)/60)/2)
times = n.arange(start_jd,end_jd,(1./24/60))
print 'Observation duration:', start_jd, end_jd

#other observing parameters
fq = .150
#dish_size_in_lambda = 5.51 #PAPER primary beam model for HERA
#dish_size_in_lambda = 7. #simple lambda/D
#sigma = 0.0727
#sigma = dish_size_in_lambda/2.35 #in uv plane
#sigma = 1/(2*n.pi*0.0727) #analytic with right Fourier convention
#sigma = 1/(2*n.pi*0.0643)

#load array information
aa = a.cal.get_aa(opts.cal,n.array([.150]))
cat = a.src.get_catalog(opts.cal,'z')
nants = len(aa)

cnt = 0
uvbins = {}
aa.set_jultime(cen_jd)
obs_lst = aa.sidereal_time()
obs_zen = a.phs.RadioFixedBody(obs_lst,aa.lat)
obs_zen.compute(aa)
for i in xrange(nants):
    print 'working on antenna %i of %i' % (i, len(aa))
    for j in xrange(nants):
        if i == j: continue
        u,v,w = aa.gen_uvw(i,j,src=obs_zen)
        uvbin = '%.1f,%.1f' % (u,v)
        if False:
            #hack for danny
            if uvbin not in ['-15.0,0.0','15.0,0.0','-15.0,-0.0','15.0,-0.0','15.0,2.0','15.0,-2.0','-15.0,2.0','-15.0,-2.0']: continue
            cnt +=1
            print cnt, i,j, uvbin
        if not uvbins.has_key(uvbin): uvbins[uvbin] = ['%i,%i' % (i,j)]
        else: uvbins[uvbin].append('%i,%i' % (i,j))

print 'There are %i baseline types' % len(uvbins.keys())
#print uvbins.keys()

dim = n.round(SIZE/dish_size_in_lambda/2)*2 - 1 # round to nearest odd
print 'dim', dim
uvplane = {}
uvsum,quadsum = n.zeros((dim,dim)), n.zeros((dim,dim))
for cnt, uvbin in enumerate(uvbins):
    print 'working on %i of %i uvbins' % (cnt+1, len(uvbins))
    uvplane = n.zeros((dim,dim))
    for t in times:
        aa.set_jultime(t)
        lst = aa.sidereal_time()
        obs_zen.compute(aa)
        bl = uvbins[uvbin][0]
        nbls = len(uvbins[uvbin])
        i, j = bl.split(',')
        i, j = int(i), int(j)
        #print i,j,bl,nbls
        u,v,w = aa.gen_uvw(i,j,src=obs_zen)
        _beam = beamgridder(sigma=sigma/dish_size_in_lambda,xcen=u/dish_size_in_lambda,ycen=v/dish_size_in_lambda,size=dim)
        #print sigma/dish_size_in_lambda, u/dish_size_in_lambda,v/dish_size_in_lambda, dim       
        uvplane += nbls*_beam
        uvsum += nbls*_beam
    quadsum += (uvplane)**2

#uv[:,:dim/2] = 0
#uv[dim/2:,dim/2] = 0

quadsum = quadsum**.5

n.savez('uvcov.npz',sum=uvsum,quadsum=quadsum)

print 'there are %i minutes of integration in the uv plane' % n.sum(uvsum)

p.imshow(quadsum,interpolation='nearest')
p.colorbar()
p.show()

if False:
    for uv in uvplane:
        p.imshow(uvplane[uv],interpolation='nearest')
        p.title(uv)
        p.colorbar()
        p.show()

