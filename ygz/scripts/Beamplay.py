import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
X,Y,Z=img.get_top(center=(0,0))
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
for i in range(nants):
  for j in range(i+1,nants):
     # for time in n.arange(2456240.2,2456241,.01):
          #  aa.set_jultime(time)
            #print aa.get_jultime()
            print 'Computing Antenna Pair', i, j
            src.compute(aa)

            ntop=n.array([X,Y,Z])
            bm1x=aa[i].bm_response(ntop,pol='x')[0]
            bm2x=aa[j].bm_response(ntop,pol='x')[0]

            bm=bm1x*n.conj(bm2x)
            bmsq=(bm1x*n.conj(bm2x))**2
            print bmsq.shape
            bmp=bmsq
            
            bmp.shape=shape0
            
            print bmp.shape
            fbm=n.fft.fft2(bmp)
            freq=n.fft.fftfreq(400,d=d)
            fbmamp=n.log(n.abs(fbm))
            
            rax=[-1,1,-1,1]
            freqax=[n.amin(freq),n.amax(freq),n.amin(freq),n.amax(freq)]

            bmp=bmp
            p.axes(rax)
            im = plt.imshow(bmp, interpolation=None)
            
                    #p.plot(time, aa.sidereal_time(),'k.')
                    #p.plot(u,v,'.',ms=2,color = (((j*17)%33)/33., (32-i)*17%33/33., (i+j)/65., 1)) 
                    #print u,v
            
            #p.plot(-u,-v,'ko')
plt.colorbar()
p.xlabel('u',size=14)
p.ylabel('v',size=14)
plt.show()
