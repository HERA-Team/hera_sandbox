#plots configuration of antenna array
import aipy as a, numpy as n, pylab as p, ephem as e
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import matplotlib.ticker as tic
aa = a.cal.get_aa('psa6240_v003',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants = 64
rad2deg = 180/n.pi
ltsec = 299792458.  #meters
X,Y,Z,I = [],[],[],[]

for i in range(nants):
    #print i
    a = aa.ants[i]
    pos = a.pos*ltsec*1.E-9       #X,Y,Z position in meters
    X.append(pos[0]); Y.append(pos[1]); Z.append(pos[2]); I.append(i)

Xt,Yt=X,Y
X,Y=Yt,Xt
fig = p.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(211)
p.scatter(X,Y)
ax.set_aspect(1)
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(212)
p.scatter(X,Y)
for x,y,i in zip(X, Y,I):
    ax.annotate('%s' %i, xy=(x,y), textcoords='data') # <--
ax.set_xlabel('East Position [m]')
ax.set_ylabel('North Position [m]')
ax.set_aspect(5)

p.grid()
p.show()
