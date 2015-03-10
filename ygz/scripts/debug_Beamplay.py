import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import Fbeam,Rbeam
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

p.figure()
aa.set_jultime(2456240.2)
time_range=n.arange(2456240.2,2456240.6,.05)
peak=[]
for i in range(nants):
  for j in range(i+1,nants):
      #for time in time_range:
       #     aa.set_jultime(time)
            #print aa.get_jultime()
            print 'Computing Antenna Pair', i, j

            ntop=n.array([X,Y,Z])
      #note all beams are the same
            bmp= Rbeam(aa[i], ntop, shape0, 'x')
            freq, fbmamp= Fbeam(bmp, d, 400)

            rax=[-1,1,-1,1]
            freqax=[freq[0],freq[len(freq)-1], freq[0],freq[len(freq)-1]]
            mid=len(freq)/2
            peak.append(fbmamp[mid,mid])
            #for dat in fbmamp:
            #    if isinstance(dat, float) != 'True':
            #        print dat
            im = plt.imshow(bmp, interpolation=None, extent=rax)
            #im = plt.imshow(fbmamp, interpolation='nearest', extent=freqax)

                    #p.plot(time, aa.sidereal_time(),'k.')
                    #p.plot(u,v,'.',ms=2,color = (((j*17)%33)/33., (32-i)*17%33/33., (i+j)/65., 1)) 
                    #print u,v
            
            #p.plot(-u,-v,'ko')
plt.colorbar()
#p.xlabel('u',size=14)
#p.ylabel('v',size=14)
mid=len(freq)/2
#plt.plot(freq[(mid-10):(mid+10)],fbmamp[mid-5,(mid-10):(mid+10)])
#plt.plot(time_range,peak)
plt.show()
