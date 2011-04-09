#mixes AIPY and CASA to grid down uv samples using central freq to compute 
#uvw
# DCJ 2011
#sick socorro suck
import numpy as n
uvsize=200
uvres =0.4 
DIM = int(uvsize/uvres)
ulim = [1.,2000.]
highchan = 400
lowchan = 300

c = 3e8 #m/s

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
print "loading data...",
ms.open(vis)
ms.selectinit()
ms.select({'uvdist':ulim})
#ms.iterinit(columns=['TIME'])
ms.selectchannel(highchan-lowchan,lowchan,1,1)
rec = ms.getdata(['axis_info'])
f = n.median(rec['axis_info']['freq_axis']['chan_freq'].squeeze())

#moredata=True
#while(moredata):
rec = ms.getdata(['u','v','w','data','flag'])
D = n.ma.array(rec['data'],mask=rec['flag'])
D = D.squeeze()

U,V,W = rec['u']*f/c,rec['v']*f/c,rec['w']*f/c
del(rec)
ms.close()
print " done"
U.shape += (1,)
V.shape += (1,)
W.shape += (1,)



print "gridding.... ",
uvs = []
for d in D:
    im = a.img.Img(uvsize, uvres, mf_order=0)
    grid_it(im,U,V,W,D,n.ones_like(D))
    print ".",
    UV = a.img.recenter(n.abs(im.uv).astype(n.float), (DIM/2,DIM/2))
    uvs.append(UV)
    del(im)
print "done"
#output the result in a casa image
#matshow(uvs[:,:,0,0])
uvs = n.array(uvs)
#output correctly to a casa image file.

