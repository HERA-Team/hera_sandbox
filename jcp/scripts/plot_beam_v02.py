#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse, pickle
from mpl_toolkits.basemap import Basemap

class DataQuery:
    def __init__(self,fig,h,pkl=None,npz=None):
        self.x = 0
        self.y = 0
        self.cid = fig.canvas.mpl_connect('button_press_event', self)
    def __call__(self,event):
        try:
        #if True:
            self.x = event.xdata - 1
            self.y = event.ydata - 1
            self.z = n.abs(n.lib.scimath.sqrt(1 - self.x**2 - self.y**2))
            crd, = h.crd2px(n.array([self.x]),n.array([self.y]),n.array([self.z]))
            print 'beam= ',h[crd]
            if pkl == None:
                print 'pix= ',crd
            else:
                pkl_file = open(pkl,'rb')
                data = pickle.load(pkl_file)
                #print data['meas'].keys()
                print 'pix= ',crd
                for k in data['meas'].keys():
                    srcs = n.load(npz)['srcnames']
                    fluxes = 10**n.load(npz)['srcfluxes']
                    if data['meas'][k].has_key(crd):
                        print k,'\tcat val= ',fluxes[n.where(srcs==k)],'\tpJys= ',data['meas'][k][crd],'\tbeam= ',(data['meas'][k][crd]/fluxes[n.where(srcs==k)]),'\twgt= ','%e' % data['wgt'][k][crd]
        except(TypeError):
            pass

o = optparse.OptionParser()
o.set_usage('plot_beam.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True,cal=True)
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('--pkl',dest='pkl',default=None,
    help="Return underlying data values from this pkl file when map is clicked.")
o.add_option('-m','--mode',dest='mode',default='linear',
    help="Scale to plot the beam in.")
opts,args = o.parse_args(sys.argv[1:])

fig = p.figure()
axis = fig.add_subplot(111)

cmap = p.get_cmap(opts.cmap)
aa = a.cal.get_aa(opts.cal,n.array([.150]))

map = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()

lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
y,x,z = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
ax,ay,az = a.coord.latlong2xyz(n.array([aa.lat,0]))
#x = x - ax
#y = y - ay
#z = z - az
try: data, indices = h[x,y,z]
except(ValueError): data = h[x,y,z]
data.shape = lats.shape

map.drawmapboundary()
map.drawmeridians(n.arange(0, 360, 30))
map.drawparallels(n.arange(0, 90, 10))

if opts.mode == 'log':
    data = n.where(data == 0.0, 0, n.log(n.abs(data)))

if opts.max is None: max = data.max()
else: max = opts.max
if opts.drng is None:
    min = data.min()
#    if min < (max - 10): min = max-10
else: min = max - opts.drng
step = (max - min) / 10
levels = n.arange(min-step, max+step, step)
print min,max
#data = data.clip(min, max)
#data = n.ma.array(data, mask=mask)
max,min=-1,-3.5
map.imshow(data, vmax=max, vmin=min, cmap=cmap)
#map.contourf(cx,cy,data,levels,linewidth=0,cmap=cmap)
p.colorbar()

pkl = opts.pkl
npz = opts.pkl.replace('pkl','npz')

query = DataQuery(fig,h,pkl=pkl,npz=npz)
#print crds.x,crds.y
#fig.canvas.mpl_disconnect(cid)

p.show()
