#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys, capo as C
#import _beamuv as beamuv
import beamuv

o = optparse.OptionParser()
o.add_option('-b','--beam',dest='beam',default=None,
    help='The beam fits file to use.')
o.add_option('-t','--text',dest='text',default=None,
    help='The text file of source fluxes to use.')
opts,args = o.parse_args(sys.argv[1:])

#_coeffs = beamuv.coeffs_from_file(opts.beam)
#beam = beamuv.BeamUV(_coeffs,.150,size=1000.,pol='y')

beam = a.map.Map(fromfits=opts.beam,interp=False)

srctimes,srcfluxes,x,y,z = C.jcp.read_srcnpz(args, verbose=True)
for src in srcfluxes: srcfluxes[src] = n.mean(n.abs(srcfluxes[src]), axis=1)
srcs = srctimes.keys()

fluxes = n.loadtxt(opts.text,usecols=[2])
srcnames = n.loadtxt(opts.text,usecols=[1],dtype='string')
fluxdict = {}
for f,k in zip(fluxes,srcnames):
    fluxdict[k] = f

for k in srcs:
    track = []
    for ix,iy,iz in zip(x[k],y[k],z[k]):
        #track.append(beam.response(n.array([ix]),n.array([iy]),n.array([iz]))**2)
        track.append(beam[ix,iy,iz])
    track = n.array(track).squeeze()
    #print track.shape

    if False:
        p.semilogy(srctimes[k],srcfluxes[k],',',label='data')
        p.semilogy(srctimes[k],track*fluxdict[k],',',label='prediction')
        p.legend()
        p.title(k)
        p.show()

    flux, hi_flux, lo_flux = fluxdict[k], fluxdict[k], fluxdict[k]

    valid = n.where(track != 0, 1, 0)
    srcfluxes[k] = srcfluxes[k].compress(valid)
    track = track.compress(valid)

    wgt = n.ones_like(track)
    var = n.sum((srcfluxes[k] - track*flux)**2)/n.sum(1./wgt)
    #print var
    wgt = 1./var
    #chi2 = n.sum(wgt*(srcfluxes[k] - track*flux)**2)/n.sum(wgt)
    chi2 = (1./(len(track)-1)) * n.sum((srcfluxes[k] - track*flux)**2/(1./wgt))
    #print chi2
    while True:
        hi_flux += .1
        #_chi2 = n.sum(wgt*(srcfluxes[k] - track*hi_flux)**2)/n.sum(wgt)
        #_chi2 = n.sum((srcfluxes[k] - track*hi_flux)**2)/n.sum(1./wgt)
        _chi2 = (1./(len(track)-1)) * n.sum((srcfluxes[k] - track*hi_flux)**2/(1./wgt))
        #if _chi2 > 2 * chi2: break
        #print _chi2-chi2
        if _chi2 > chi2 + 1.: break       
    while True:
        lo_flux -= .1
        #_chi2 = n.sum(wgt*(srcfluxes[k] - track*lo_flux)**2)/n.sum(wgt)
        #_chi2 = n.sum((srcfluxes[k] - track*lo_flux)**2)/n.sum(1./wgt)
        _chi2 = (1./(len(track)-1)) * n.sum((srcfluxes[k] - track*lo_flux)**2/(1./wgt))
        #if _chi2 > 2 * chi2: break
        #if _chi2 < chi2: print _chi2-chi2
        if _chi2 > chi2 + 1.: break      
    print k, flux, '+', hi_flux-flux, '-', flux-lo_flux
