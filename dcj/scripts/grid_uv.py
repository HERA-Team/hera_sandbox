#mixes AIPY and CASA to grid down uv samples using central freq to compute 
#uvw
# DCJ 2011
#sick socorro suck
import numpy as n
from cosmo_units import *
uvsize=200
uvres =4 
umin = 45
umax = 200
nu = 6
DIM = int(uvsize/uvres)
ulim = [1.,2000.]
highchan = 450
lowchan = 300

c = 3e8 #m/s
from numpy import array,int32
csystemplate = {'linear0': {'axes': array(['UU', 'VV'], 
      dtype='|S3'),
             'cdelt': array([-2.29183118,  2.29183118]),
             'crpix': array([ 250.,  250.]),
             'crval': array([ 0.,  0.]),
             'pc': array([[ 1.,  0.],
       [ 0.,  1.]]),
             'units': array(['lambda', 'lambda'], 
      dtype='|S7')},
 'linear2': {'axes': array(['Time'], 
      dtype='|S5'),
             'cdelt': array([  2.04834084e-07]),
             'crpix': array([ 50.]),
             'crval': array([ 0.]),
             'pc': array([[ 1.]]),
             'units': array(['s'], 
      dtype='|S2')},
 'obsdate': {'m0': {'unit': 'd', 'value': 55455.173560246825},
             'refer': 'UTC',
             'type': 'epoch'},
 'observer': '',
 'pixelmap0': array([0, 1], dtype=int32),
 'pixelmap1': array([2], dtype=int32),
 'pixelmap2': array([3], dtype=int32),
 'pixelreplace0': array([ 0.,  0.]),
 'pixelreplace1': array([ 0.]),
 'pixelreplace2': array([ 0.]),
 'pointingcenter': {'initial': False,
                    'value': array([ 1.3962    , -0.53619133])},
 'stokes1': {'axes': array(['Stokes'], 
      dtype='|S7'),
             'cdelt': array([ 1.]),
             'crpix': array([ 0.]),
             'crval': array([ 1.]),
             'pc': array([[ 1.]]),
             'stokes': array(['I'], 
      dtype='|S2')},
 'telescope': 'PAPER',
 'worldmap0': array([0, 1], dtype=int32),
 'worldmap1': array([2], dtype=int32),
 'worldmap2': array([3], dtype=int32),
 'worldreplace0': array([ 0.,  0.]),
 'worldreplace1': array([ 1.]),
 'worldreplace2': array([ 0.])
}



def grid_it(im,us,vs,ws,ds,wgts):
    #print 'Gridding %d integrations' % n_ints
    sys.stdout.write('|'); sys.stdout.flush()
    if len(ds) == 0: raise ValueError('No data to use.')
    ds,wgts = n.concatenate(ds), n.concatenate(wgts).flatten()
    us,vs,ws = n.concatenate(us), n.concatenate(vs), n.concatenate(ws)
    # Grid data into UV matrix
    (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
    im.put((us,vs,ws), ds, wgts)
    #im.put((us,vs,ws), ds, wgts, invker2=bm_im)
def plop(X,x):
    return n.argwhere(n.min(n.abs(X-x))==n.abs(X-x)).squeeze()
def length(X): return n.sqrt(n.dot(X,X))
def MJDs2JD(t):
    return t/86400 + 2400000.5
def JD2MJDs(t):
    return (t-2400000.5)*86400

flush = sys.stdout.flush
print "loading data...",
flush()
ms.open(vis)
ms.selectinit()
ms.select({'uvdist':ulim})
#ms.iterinit(columns=['TIME'])
ms.selectchannel(highchan-lowchan,lowchan,1,1)
rec = ms.getdata(['axis_info'])
f = n.median(rec['axis_info']['freq_axis']['chan_freq'].squeeze())
df = rec['axis_info']['freq_axis']['resolution'].squeeze()[0]
#moredata=True
#while(moredata):
rec = ms.getdata(['u','v','w','data','flag','antenna1','antenna2','time'])
t = n.median(rec['time'])
I,J = rec['antenna1'],rec['antenna2']
D = n.ma.array(rec['data'],mask=rec['flag'])
D = D.squeeze()

U,V,W = rec['u']*f/c,rec['v']*f/c,rec['w']*f/c
del(rec)
ms.close()
tb.open(vis+'/FIELD')
direction = tb.getcol('PHASE_DIR').squeeze()
tb.close()



print " done"
flush()

print "match uvw with ijs",;flush()
Nant = n.max(n.vstack((I,J)))
UVWs = n.zeros((Nant,Nant,3))
bls = {}
#for i,j in n.indices(UVWs.shape):
for (i,j,u,v,w) in zip(I,J,U,V,W):
    bl = '%d&&%d'%(i,j)
    if not bls.has_key(bl): bls[bl] = []
    bls[bl].append((u,v,w))
for bl in bls:
    BL = n.array(bls[bl])
    bls[bl] = (n.median(BL[:,0]),n.median(BL[:,1]),n.median(BL[:,2]))
print ". done";flush()
#TODO grid the single set of uvws, use get_indices()
#TODO Then grid by the incoherent uvs radially.  
ci = [] #coherent indices (into Im a uv grid)
ici = [] #incoherent indices (into Ps, a list of radial us)

Im = a.img.Img(uvsize, uvres, mf_order=0)
Ps = n.logspace(n.log10(umin),n.log10(umax),num=(nu+2))  #we'll eventually throw out 0 & -1
print "badly written, but exact, gridding",;flush()
for i,j in zip(I,J):
    bl = '%d&&%d'%(i,j)
    u,v = bls[bl][0],bls[bl][1]
    uv = n.array([u,v])
    if length(uv)<Ps.min() or length(uv)>Ps.max(): 
        ci.append((-1,-1))
        ici.append(-1)
    else:
        ci.append(Im.get_indices(u,v)
)
        ici.append(plop(Ps,length(uv)))
uvs = n.zeros(Im.uv.shape+(D.shape[0],))    #uvf cube
uvi = n.zeros(Im.uv.shape)-1                  #map to radial bins
uvin = n.zeros_like(uvs).astype(n.int)
for l,(ui,vi) in enumerate(ci):
    if ui<0:continue
    uvs[ui,vi,:] += D[:,l]
    uvi[ui,vi] = ici[l]
    uvin[ui,vi,:] += 1
print ".. done";flush()
uvs[uvin>1] /= uvin[uvin>1]
#STOP! Save the uvgrid and exit.
print "FFT",flush()
nchan = uvs.shape[2]
uveta = n.abs(n.fft.fft(uvs,axis=2))[:,:,:nchan/2]
uveta = a.img.recenter(uveta,(uvs.shape[0]/2,uvs.shape[1]/2,0))
dly = n.fft.fftfreq(uveta.shape[2],df)[:nchan/2]
csupdate = {'linear0':{'cdelt':n.array([-uvres, uvres]).astype(n.float),
                        'crpix':n.array((uveta.shape[0]/2,uveta.shape[1]/2)).astype(n.float)},
            'obsdate':{'m0':{'unit':'d','value':t/86400},
                        'refer':'UTC','type':'epoch'},
            'pointingcenter':{'initial':False,
                                'value':direction
                       },
            'linear2':{'cdelt':array([n.diff(dly)[0]]),
                        'crpix':array([0.]),
                        'crval':array([dly[0]])}    
                }
csysrecord = csystemplate
for k in csupdate:
    csysrecord[k].update(csupdate[k])

#csys = cs.fromrecord(csysrecord)
print "writing ",vis+'.uveta';flush()
uveta.shape = (uveta.shape[0],uveta.shape[1],1,uveta.shape[2])
ia.fromarray(outfile=vis+'.uveta',pixels=uveta,csys=csysrecord,overwrite=True)
