#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

calfile = 'psa898_v002' 
na = 3#2 # number of antennas in array
freqs = n.arange(.1,.2,.01) # frequencies to loop over

aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array

im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
valid = n.logical_not(tx.mask)
tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()

#beam response for an antenna pointing at (tx,ty,tz) with a polarization in x direction
#amp = A(theta) in notes
amp = aa[0].bm_response((tx,ty,tz),pol='x')**2 

#coefficient array: rows baselines; cols frequencies
coeffs = n.zeros([(na*na-na)/2,len(freqs)])

baselines = n.zeros([(na*na-na)/2,3])
# loop through all baselines
ll=0
for ii in n.arange(na):
    for jj in n.arange(ii):
        bx,by,bz = aa.get_baseline(ii,jj,'z') 
        #the baseline for antennas ii and jj 
        print ii,jj,[bx,by,bz]
        baselines[ll] = [bx,by,bz]

        kk=0
        for fq in freqs: #loop over frequencies
            phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
            phs.shape = amp.shape = im.uv.shape #?????
            amp = n.where(valid, amp, 0)
            phs = n.where(valid, phs, 0)

            dc_response = n.sum(amp*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
            jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
            print '\t',fq, dc_response, jy_response
            coeffs[ll,kk] = dc_response
            kk+=1
        ll+=1
n.savetxt('{0}_baselines.txt'.format(calfile),baselines)
n.savetxt('{0}_frequencies.txt'.format(calfile),freqs)
n.savetxt('{0}_coeffs.txt'.format(calfile),coeffs)

