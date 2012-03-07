#! /usr/bin/env python
import numpy as n, pylab as pl, aipy as a, optparse, sys

o = optparse.OptionParser()
o.set_usage('pipeline.py [uvfile(s)]')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True,pol=True,ant=True)
o.add_option('--size', dest='size', type='int', default=300,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
opts, args = o.parse_args(sys.argv[1:])

def dtransform(vis):
    #the delay transform on one miriad visibility array
    _vis = n.fft.fft(vis)
    return _vis

#def _vis2pk():
    #turns delay transformed visibilties to power spectrum

im = a.img.Img(opts.size,opts.res)
DIM = int(opts.size/opts.res)


us,vs,ws,ds,wgts = [],[],[],[],[]
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'],uv['nchan'])

#phase to the first integration of the first file
(crd,t,(i,j)),d = uv.read()
aa.set_jultime(t)
s = a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, name='zen')
src = a.fit.SrcCatalog([s])

for file in args:
    uv = a.miriad.UV(file)
    a.scripting.uv_selector(uv,opts.ant,opts.pol)
    for p,d,f in uv.all(raw=True):
        crd,t,(i,j) = p
        aa.set_jultime(t)
        src.compute(aa)
        #d /= aa.passband(i,j)
        
        wgt = n.ones(d.shape,dtype=n.float)
        valid = n.logical_not(f)
        d = d.compress(valid)
        if len(d) == 0: continue
        u,v,w = aa.gen_uvw(i,j,src=s)
        us.append(u.compress(valid))
        vs.append(v.compress(valid))
        ws.append(w.compress(valid))
        ds.append(d)
        wgts.append(wgt.compress(valid))

#grid for plotting
ds,wgts = n.concatenate(ds), n.concatenate(wgts).flatten()
us,vs,ws = n.concatenate(us), n.concatenate(vs), n.concatenate(ws)
(us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
#im.put((us,vs,ws), ds, wgts)
im.put((us,vs,ws),ds)
uvs = a.img.recenter(n.abs(im.uv).astype(n.float), (DIM/2,DIM/2))
uvs = n.where(uvs > 0., n.log10(uvs),-20)

pl.imshow(uvs)
pl.colorbar()
pl.show()
