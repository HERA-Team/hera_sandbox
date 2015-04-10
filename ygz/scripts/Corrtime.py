import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import interpolate
import export_beam

sz = 200
d = 1./sz
img = a.img.Img(200,res=0.5)
#400 by 400 image, i.e. 200 boxes, with 4 pixels per box
#this will give freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
aa = a.cal.get_aa('psa6622_v001',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

ant1 = 114
ant2 = 115
ant3 = 24
ant4 = 113

aa.set_jultime(2456240.3)
src.compute(aa)
#u0,v0,w0 = aa.gen_uvw(ant1,ant2,src=src)
#u0,v0 = u0.flatten(), v0.flatten()
U = []
V = []
times = n.arange(2456240.2,2456240.3,.001)
for clock in times:
    aa.set_jultime(clock)
    src.compute(aa)
    u0,v0,w0 = aa.gen_uvw(ant1,ant2,src=src)
    u0,v0 = u0[0][0], v0[0][0]
    if u0<0: u0,v0 = -u0,-v0
    u1,v1,w1 = aa.gen_uvw(ant3,ant4,src=src)
    u1,v1 = u1[0][0], v1[0][0]
    if u1<0: u1,v1 = -u1,-v1
    u,v = u0-u1,v0-v1
    U.append(u)
    V.append(v)


U,V = n.array(U), n.array(V)
U,V = U.flatten(), V.flatten()
ntop = n.array([X,Y,Z])
bmp = export_beam.beam_real(aa[ant1], ntop, shape0, 'x')
freq, fbmamp = export_beam.beam_fourier(bmp, d, 400)
val = export_beam.get_overlap(freq,fbmamp,U,V)

fig = p.figure()
ax1 = fig.add_subplot(121)
p.plot(times,val)
ax1.set_xlabel('jultime',size=12)
ax1.set_ylabel('FT(A^2)',size=12)
#p.plot(TIME,U)
#p.plot(TIME,V)
ax2 = fig.add_subplot(122)
p.plot(U,V)
ax2.set_xlabel('U',size=12)
ax2.set_ylabel('v',size=12)
plt.show()
