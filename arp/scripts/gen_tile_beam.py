#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--d', type='float', help='Element spacing in wavelengths.')
o.add_option('--nelem', type='int', help='Number of elements on a tile side.')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))

im = a.img.Img(size=1000, res=.5)
#im = a.img.Img(size=10, res=.5)

NELEM = opts.nelem
D = opts.d
#N = 2.0 * a.const.c / .150e9 
N = D * a.const.c / .150e9 

antpos = {}
cnt = 0
grid = n.arange(-float(NELEM)/2,float(NELEM)/2) + 0.5
#print N * grid
for i in N * grid:
    for j in N * grid:
        antpos['%d'%cnt] = {'top_x': i, 'top_y': j, 'top_z': 0.}
        cnt += 1
aa.set_params(antpos)

u,v = [],[]

for i in range(cnt):
    for j in range(i, cnt):
        _u,_v,_w = aa.gen_uvw(i, j)
        u.append(_u); v.append(_v)
        if i != j:
            u.append(-_u); v.append(-_v)

u,v = n.array(u).squeeze(), n.array(v).squeeze()
d,w = n.ones_like(u), n.zeros_like(u)

#for _u,_v,_d in zip(u,v,d):
#    print _u,_v,_d
im.put((u,v,w), d)

img = im.image()

x,y,z = im.get_top()
x,y,z = x.flatten(), y.flatten(), z.flatten()
img = img.flatten()
resp = aa[0].bm_response((x,y,z))[0]**2
resp = n.where(z <= 0, 0, resp)
beam = resp * n.abs(img)

h = a.map.Map(nside=64, interp=False)
h.add((x,y,z), n.ones_like(beam), n.sqrt(beam)) # writing a voltage beam
h.reset_wgt()
filename = 'beam_%dx%d_d%3.1f.hmap' % (NELEM,NELEM,D)
print 'Writing', filename
h.map.to_fits(filename)

if False:
    beam.shape = im.shape
    p.subplot(121)
    p.imshow(a.img.recenter(beam, n.array(im.shape)/2))
    p.colorbar(shrink=.5)
    #p.imshow(im.image(n.array(im.shape)/2))
    p.subplot(122)
    p.imshow(a.img.recenter(n.abs(im.uv), n.array(im.shape)/2))
    p.colorbar(shrink=.5)
    p.show()
