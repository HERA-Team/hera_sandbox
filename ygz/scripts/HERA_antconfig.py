#plots configuration of antenna array
import aipy as a, numpy as n, pylab as p, ephem as e
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import matplotlib.ticker as tic
import seaborn as sns

FILE = "../calfiles/HERA_antconfig/antenna_positions_350.dat"

antpos = np.loadtxt(FILE)
nants = antpos.shape[0]
X,Y,Z = antpos[:,0], antpos[:,1], antpos[:,2]
X -= n.average(X)
Y -= n.average(Y)
I = n.arange(nants)

sns.set(style="darkgrid")
sns.set_context("poster")
fig = p.figure()
# #ax = fig.gca(projection='3d')
# ax = fig.add_subplot(211)
# p.scatter(X,Y)
# ax.set_aspect(1)
# setp( ax.get_xticklabels(), visible=False)
# setp( ax.get_yticklabels(), visible=False)

# ax = fig.add_subplot(212)
ax = fig.add_subplot(111)
g = 320
Xg, Yg, Ig = X[:g], Y[:g], I[:g]
p.scatter(Xg,Yg)
for x,y,i in zip(Xg, Yg,Ig):
    ax.annotate('%s' %i, xy=(x,y), textcoords='data', fontsize=12) # <--
ax.set_xlabel('East Position [m]')
ax.set_ylabel('North Position [m]')
ax.set_aspect(1)
#ax.grid

p.tight_layout()
p.show()
