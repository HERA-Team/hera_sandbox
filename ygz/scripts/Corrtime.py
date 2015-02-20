import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import interpolate
sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
#400 by 400 image, i.e. 200 boxes, with 4 pixels per box
#this will give freq space kmax=100, dk=0.5
X,Y,Z=img.get_top(center=(200,200))
shape0=X.shape
X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()

aa=a.cal.get_aa('psa898_v003',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants=32
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

ant1=9
ant2=10
ant3=15
ant4=16

aa.set_jultime(2456240.3)
src.compute(aa)
u0,v0,w0 = aa.gen_uvw(ant1,ant2,src=src)
u0,v0 = u0.flatten(), v0.flatten()

print 'u0v0', u0, v0
p.figure()

U=[]
V=[]
TIME=n.arange(2456240.2,2456240.6,.01)
for time in TIME:
    aa.set_jultime(time)
    #print aa.get_jultime()
    src.compute(aa)


    u1,v1,w1 = aa.gen_uvw(ant3,ant4,src=src)
    u1,v1 = u1.flatten(), v1.flatten()

    u=u0-u1
    v=v0-v1
    U.append(u)
    V.append(v)


U, V=n.array(U), n.array(V)
U, V=U.flatten(), V.flatten()

ntop=n.array([X,Y,Z])
bm1x=aa[ant1].bm_response(ntop,pol='x')[0]
bm2x=aa[ant2].bm_response(ntop,pol='x')[0]
bm3x=aa[ant3].bm_response(ntop,pol='x')[0]
bm4x=aa[ant4].bm_response(ntop,pol='x')[0]
bm=bm1x*n.conj(bm2x)
bmsq=(bm1x*n.conj(bm2x))*(bm3x*n.conj(bm4x))
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
f = interpolate.interp2d(freq, freq, fbmamp, kind='cubic')


val=n.diagonal(f(U,V))

#p.plot(TIME,val)
#p.plot(TIME,U)
#p.plot(TIME,V)
p.plot(U,V)
p.xlabel('jultime',size=12)
p.ylabel('FT(A^2)',size=12)
plt.show()
