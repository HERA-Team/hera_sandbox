__author__ = 'yunfanzhang'
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

aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants=32
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

ant1=0
ant2=12
antdict={}

dt=0.003
TIME=n.arange(2456240.2,2456240.3, dt)

distlim=0.01

for time in TIME:
    aa.set_jultime(time)
    src.compute(aa)
    u0,v0,w0 = aa.gen_uvw(ant1,ant2,src=src)

    u0,v0 = u0[0,0], v0[0,0]

    antdict[u0,v0]=[[ant1,ant2,time]]

for i in range(nants):
    print i
    for j in range(i+1,nants):
        for time in TIME:
            aa.set_jultime(time)
            src.compute(aa)
            u1,v1,w1 = aa.gen_uvw(i,j,src=src)
            u1,v1 = u1.flatten()[0], v1.flatten()[0]
            for key in antdict.keys():
                u=key[0]-u1
                v=key[1]-v1
                #print u,v
                delt=n.abs(time-antdict[key][0][2])
                for entry in antdict[key]:
                    if n.abs(entry[2]-time)<delt:
                        delt=n.abs(entry[2]-time)

                if n.sqrt(u**2+v**2)<distlim and delt>1.5*dt:
                    antdict[key].append([i,j,time])
                    #print key
                else:
                    antdict[u1,v1]=[[i,j,time]]


for key in antdict.keys():
    if len(antdict[key])<2:
        del antdict[key]
print antdict

