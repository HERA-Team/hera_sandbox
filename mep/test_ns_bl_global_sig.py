#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

aa = a.cal.get_aa('psa898_v002', n.array([.150])) #get antenna array

im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
valid = n.logical_not(tx.mask)
tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()

amp = aa[0].bm_response((tx,ty,tz),pol='x')**2 
#this is the beam response for an antenna pointing at (tx,ty,tz) 
#with a polarization in x direction
#amp = A(theta) in notes

bx,by,bz = aa.get_baseline(0,1,'z') 
#the baseline for antennas 0 and 1 in the array
print [bx,by,bz]

b2x,b2y,b2z = aa.get_baseline(1,28,'z') 
#the baseline for antennas 1 and 28 in the array
print [b2x,b2y,b2z]

for fq in n.arange(.1,.2,.01): #loop over frequencies
    phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
    phs.shape = amp.shape = im.uv.shape #?????
    amp = n.where(valid, amp, 0)
    phs = n.where(valid, phs, 0)

    dc_response = n.sum(amp*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
    jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
    print fq, dc_response, jy_response

    phs2 = n.exp(-2j*n.pi*fq * (b2x*tx+b2y*ty+b2z*tz)) #fringe pattern
    phs2.shape = amp.shape = im.uv.shape #?????
    amp = n.where(valid, amp, 0)
    phs2 = n.where(valid, phs2, 0)

    dc_response2 = n.sum(amp*phs2)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
    jy_response2 = n.real(dc_response2 * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
    print fq, dc_response2, jy_response2

    if True:
        p.subplot(231); p.imshow(amp)
        p.subplot(232); p.imshow(amp*phs.real,vmax=1,vmin=-1)
        p.subplot(233); p.imshow(phs.real)

        p.subplot(234); p.imshow(amp)
        p.subplot(235); p.imshow(amp*phs2.real,vmax=1,vmin=-1)
        p.subplot(236); p.imshow(phs2.real)

        p.show()

