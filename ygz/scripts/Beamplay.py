import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import Fbeam
sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
X,Y,Z=img.get_top(center=(200,200))
shape0=X.shape
X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()

aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
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
            bmp= Fbeam.Rbeam(aa[i], ntop, shape0, 'x')
            freq, fbmamp= Fbeam.Fbeam(bmp, d, 400)

            rax=[-1,1,-1,1]
            freqax=[freq[0],freq[len(freq)-1], freq[0],freq[len(freq)-1]]
            #mid=len(freq)/2
            #peak.append(fbmamp[mid,mid])
            #im = plt.imshow(bmp, interpolation=None, extent=rax)
            im = plt.imshow(fbmamp, interpolation='nearest', extent=freqax)

plt.colorbar()
#p.xlabel('u',size=14)
#p.ylabel('v',size=14)
plt.show()
