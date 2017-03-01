import aipy as a, numpy as np, pylab as plt, sys, os, ephem, optparse
from mpl_toolkits.basemap import Basemap

DATADIR = '../calfiles/'
XFILE = DATADIR+'GX4Y2H_4900_150.hmap'
YFILE = DATADIR+'GY4Y2H_4900_150.hmap'
map = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
Xh = a.map.Map(fromfits=XFILE)
Yh = a.map.Map(fromfits=YFILE)
for H in [Xh, Yh]:
	H.set_interpol(True)
	#H = np.where(H>0, H, 0)




# h = a.healpix.HealpixMap(nside=64)
# h.set_interpol(True)
# tx,ty,tz = h.px2crd(np.arange(h.map.size), ncrd=3)

img = a.img.Img(200,res=0.5)
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()

ant = a.cal.get_aa('psa6622_v001',np.array([.15]))[0]
PAPERBmI = ant.bm_response([X,Y,Z],pol='I').reshape((400,400))
# Bmx = np.where(Xh[X,Y,Z].real>0, Xh[X,Y,Z], 0)
# Bmx = np.where(Xh[X,Y,Z].real>0, Xh[X,Y,Z], 0)

bmI = np.sqrt((Xh[X,Y,Z].conj()*Xh[X,Y,Z] + Yh[X,Y,Z].conj()*Yh[X,Y,Z])*0.5)
#bmI = np.sqrt((Xh[X,Y,Z].conj()*Xh[X,Y,Z]*2)*0.5)
bmI /= np.amax(bmI)
bmI = bmI.reshape((400,400))

f = plt.figure()
ax = f.add_subplot(121)
map.drawmapboundary()
map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(0, 90, 10))
#import IPython; IPython.embed()

map.imshow(np.abs(bmI*bmI))
ax2 = f.add_subplot(122)
map.drawmapboundary()
map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(0, 90, 10))
#import IPython; IPython.embed()
im = map.imshow(np.abs(PAPERBmI*PAPERBmI))
f.colorbar(im, ax=[ax,ax2], shrink=0.5)
plt.show()