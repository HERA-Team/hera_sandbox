#casapy script to fool with gridded uveta power spectra
from cosmo_units import *
def length(X):
    return n.sqrt(n.dot(X,X))
def plop(X,x):
    return n.argwhere(n.min(n.abs(X-x))==n.abs(X-x)).squeeze()
def issubplotedge(rows,cols,index):
    #input same as subplot
    #return whether the subplot is on the left and or bottom edge
    # (isleft,isbottom)
    inds = n.arange(1,cols*rows+1,1)
    inds.shape = (rows,cols)
    left = not n.argwhere(inds==index).squeeze()[1]
    bottom = not (n.argwhere(inds==index).squeeze()[0]+1)%rows
#    print inds,n.argwhere(inds==index)
    return left,bottom
urmin = 40
urmax = 200
nur   = 4
band = 3 #hack to choose a channel for individual plotting

Ur = n.logspace(n.log10(urmin),n.log10(urmax),num=nur)
Urgrid = [[] for i in range(nur)]
Urgrid_dly = [[] for i in range(nur)]
psmax = psmin = 0
pl.figure(18)
pl.clf()
for file in files:
    #load file
    ia.open(file)
    uveta = ia.getchunk()
    csys = ia.coordsys()
    t = csys.epoch()['m0']['value']
#    z = float(csys.observer())
    f = csys.toworld([0,0,band,0])['numeric'][2]
    z = f212z(f)
    print f,z
    ia.close()
    #find the delay axis
    dly = [csys.toworld([0,0,band,i])['numeric'][3] for i in\
        range(uveta.shape[3])]

    UVpxs = n.argwhere(uveta[:,:,band,0]).astype(float)
    UVs = []
    for UVpx in UVpxs:
        UVworld = csys.toworld(list(UVpx)+[0,0])['numeric'][:2]
        UVs.append(UVworld)
        if Ur.min()>length(UVworld) or Ur.max()<length(UVworld):
            continue
        print UVpx,UVworld,length(UVworld),plop(Ur,length(UVworld)),
        Urgrid[plop(Ur,length(UVworld))].append(
                uveta[UVpx[0],UVpx[1],band,:]) 
        print len(Urgrid[plop(Ur,length(UVworld))])
        Urgrid_dly[plop(Ur,length(UVworld))].append(dly)
    if uveta.max()>psmax:psmax=uveta.max()
    if uveta.min()<psmin:psmin=uveta.min()

for i in n.arange(nur-1,-1,-1):
    if len(Urgrid[i])==0: del(Urgrid[i])
#do some nice plotting

rows = ceil(n.sqrt(len(Urgrid)))
cols = ceil(len(Urgrid)/rows)
nplots = rows*cols
pl.suptitle(', '.join(files)+'\n z = %4.1f'%z)
for ri,ur in enumerate(Urgrid):
    pl.subplot(rows,cols,ri+1)
    for nps,Ps in enumerate(ur):
        kparr = n.array([eta2kparr(d,z) for d in Urgrid_dly[ri][nps]])
        pl.loglog(kparr,Ps,'k',alpha=0.2)
    pl.loglog(kparr,n.mean(ur,axis=0),label=str(int(Ur[ri]))+'$\lambda$')
    pl.legend()
    pl.ylim([psmin,psmax])
    isleft,isbottom = issubplotedge(rows,cols,ri+1)
    if not isleft: pl.yticks([])
    if not isbottom: pl.xticks([])
pl.subplots_adjust(wspace=0,hspace=0)



#get u,v,ur coordinates for each 
    #sort power spectra into radial bins
    #plot each uv cell in ur subplots (we're looking for outliers)
