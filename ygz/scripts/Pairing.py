__author__ = 'yunfanzhang'
import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import interpolate
import Fpair

sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
#400 by 400 image, i.e. 200 boxes, with 4 pixels per box
#this will give freq space kmax=100, dk=0.5
X,Y,Z=img.get_top(center=(200,200))
shape0=X.shape
X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()

aa    = a.cal.get_aa('psa6622_v001',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
src   = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

nants= 128
dt   = 0.001
TIME = n.arange(2456240.2,2456240.3, dt)
dist = 0.1#size of cells to store in dictionary.

d = Fpair.Fpair_coarse(aa,nants, src,TIME,dt,dist)

for key in d.keys():
    print d[key]

fig = p.figure()
ax = fig.add_subplot(111)
U=[]
V=[]
for key in d.keys():
    for item in d[key]:
        U.append(item[2][0])
        V.append(item[2][1])
p.plot(U,V,'.',ms=5)

p.grid()
#p.xlim(-200,200)
#p.ylim(-200,200)
p.xlabel('u',size=14)
p.ylabel('v',size=14)
p.show()