import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
X,Y,Z=img.get_top(center=(200,200))
shape0=X.shape
X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()

aa=a.cal.get_aa('psa898_v003',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants=2
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

p.figure()
aa.set_jultime(2456240.2)
time_range=n.arange(2456240.2,2456240.6,.05)
peak=[]
for i in range(nants):
  for j in range(i+1,nants):
      for time in time_range:
            aa.set_jultime(time)
            #print aa.get_jultime()
            print 'Computing Antenna Pair', i, j
            src.compute(aa)

            ntop=n.array([X,Y,Z])
            bm1x=aa[i].bm_response(ntop,pol='x')[0]
            bm2x=aa[j].bm_response(ntop,pol='x')[0]

            bm=bm1x*n.conj(bm2x)
            bmsq=(bm1x*n.conj(bm2x))**2
            print bmsq.shape

            #Tranform the square of the beams
            bmp=bmsq
            bmp.shape=shape0
            
            print bmp.shape
            fbm=n.fft.fft2(bmp)
            frequv=n.fft.fftfreq(400,d=d)
            freqk=frequv*2*n.pi
            fbmamp=n.abs(fbm)
            #fbmamp=n.abs(fbm)

            numf=40

            freq=frequv
            fbmamp=n.fft.fftshift(fbmamp)
            freq=n.fft.fftshift(freq)
            #fbmamp=fbmamp[(400-numf):, :numf]
            #freq=frequv[:numf]

            rax=[-1,1,-1,1]
            freqax=[freq[0],freq[len(freq)-1], freq[0],freq[len(freq)-1]]
            mid=len(freq)/2
            peak.append(fbmamp[mid,mid])

            #im = plt.imshow(bmp, interpolation=None, extent=rax)
            #im = plt.imshow(fbmamp, interpolation=None, extent=freqax)

            
                    #p.plot(time, aa.sidereal_time(),'k.')
                    #p.plot(u,v,'.',ms=2,color = (((j*17)%33)/33., (32-i)*17%33/33., (i+j)/65., 1)) 
                    #print u,v
            
            #p.plot(-u,-v,'ko')
#plt.colorbar()
#p.xlabel('u',size=14)
#p.ylabel('v',size=14)
mid=len(freq)/2
#plt.plot(freq[(mid-10):(mid+10)],fbmamp[mid-5,(mid-10):(mid+10)])
plt.plot(time_range,peak)
plt.show()
