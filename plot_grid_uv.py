#casapy script to fool with gridded uveta power spectra
from cosmo_units import *
from shutil import move as mv
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
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
def rebin(newbins,oldbins,A):
    #returns A rebinned onto bins
    O = n.zeros_like(newbins)
    NO = n.zeros_like(newbins)
    for i,b in enumerate(oldbins):
        O[plop(newbins,b)] += A[i]
        NO[plop(newbins,b)] += 1
    O[NO>0] /= NO[NO>0]
    return O
def MJDs2JD(t):
    return t/86400 + 2400000.5
def MJD2JD(t):
    return t + 2400000.5
flush = sys.stdout.flush
urmin = 20
urmax = 150
nur   = 6
band = 3 #hack to choose a channel for individual plotting
nkparr = 10 #the number of bins to use to log-rebin in kparr direction
dologk = False
calfile = 'psa455_v004_gc'

Ur = n.logspace(n.log10(urmin),n.log10(urmax),num=nur)
Urgrid = [[] for i in range(nur)]
Urgrid_dly = [[] for i in range(nur)]
Urspec = None
Urspecn = None
psmax = psmin = 0
aa = a.cal.get_aa(calfile,n.array([0.15]))
for file in files:
    #load file
    ia.open(file)
    uveta = ia.getchunk()
    if dologk:
        Nk = nkparr
    else:
        Nk = uveta.shape[3]
    csys = ia.coordsys()
    ia.close()

    t = csys.epoch()['m0']['value']
    aa.set_jultime(MJD2JD(t))
    lst = str(aa.sidereal_time())
    fs = [csys.toworld([0,0,i,0])['numeric'][2] for i in range(uveta.shape[2])]

    if Urspec is None: 
        #assumes that all files have the same number of bands.
        Urspec = n.zeros((nur,len(fs),Nk)) 
        Urspecn = n.zeros((nur,len(fs),Nk))
        UrspecK = n.zeros((nur,len(fs),Nk))
        Urspecz = n.zeros((nur,len(fs),Nk))


    for band,f in enumerate(fs):
    #    f = csys.toworld([0,0,band,0])['numeric'][2]
        z = f212z(f)
        print f,z;flush()
        #find the delay axis
        dly = n.array([csys.toworld([0,0,band,i])['numeric'][3] for i in\
            range(uveta.shape[3])])
        #build kparr,z coordinate arrays for display purposes
        kparr = n.array([eta2kparr(d,z) for d in dly])
        logk  = n.logspace(n.log10(kparr[kparr>0].min()),
                            n.log10(kparr.max()),
                            num=nkparr)
        for i,UK in enumerate(UrspecK):
            if dologk:
                UK[band,:] = logk
            else:
                UK[band,:] = n.sqrt(kparr**2 + u2kperp(Ur[i],z)**2)
        Urspecz[:,band,:] = z
        #grid by radial U and redshift
        UVpxs = n.argwhere(uveta[:,:,band,0]).astype(float)
        UVs = []
        for UVpx in UVpxs:
            UVworld = csys.toworld(list(UVpx)+[0,0])['numeric'][:2]
            UVs.append(UVworld)
            if Ur.min()>length(UVworld) or Ur.max()<length(UVworld):
                continue
#            print UVpx,UVworld,length(UVworld),plop(Ur,length(UVworld))
            nchan = uveta.shape[3]
#            PS = (uveta[UVpx[0],UVpx[1],band,-nchan/2:] + \
#                   uveta[UVpx[0],UVpx[1],band,:nchan/2][::-1])/2
            PS = uveta[UVpx[0],UVpx[1],band,:]
            if dologk:
                ps = rebin(logk,kparr,PS)
            else:
                ps = PS
            Urspec[plop(Ur,length(UVworld)),band,:] += ps
            Urspecn[plop(Ur,length(UVworld)),band,:] += 1

        if uveta.max()>psmax:psmax=uveta.max()
        if uveta.min()<psmin:psmin=uveta.min()
    Urspec[Urspecn>0] /= Urspecn[Urspecn>0]
#    Urspec = Urspec.filled(0)
#    print type(Urspec)
    pl.figure(24)
    pl.clf()
    rows = ceil(n.sqrt(len(Ur)))
    cols = ceil(len(Ur)/rows)
    nplots = rows*cols
#    pl.suptitle(', '.join(files))
    pl.figtext(0.05,0.95,'LST:'+lst + '\n'+str(aa.epoch),family='sans-serif')
#    Urspec = n.ma.masked_where(Urspec==0,Urspec)
    cax = pl.axes([0.925,0.025,0.025,0.9])
    for i,Urs in enumerate(Urspec):
        ax = pl.subplot(rows,cols,i+1)
        if dologk:
            pl.pcolor(n.log10(UrspecK[i,:,:]),Urspecz[i,:,:],Urs,vmin=psmin,vmax=psmax)
        else:
    #        pl.pcolor(n.log10(UrspecK[i,:,:].clip(1e-2,1e9)),Urspecz[i,:,:],Ur)
            pl.pcolor(n.log10(UrspecK[i,:,:].clip(1e-2,1e9)),
                Urspecz[i,:,:],n.log10(Urs),
                vmin=n.log10(Urspec[Urspec>0].min()),vmax=n.log10(Urspec.max()))
            pl.vlines(n.log10(0.125*1.5**(n.log2(Ur[i])-4)),
                Urspecz.min(),Urspecz.max())
            print 0.125*1.5**(n.log2(Ur[i])-4)
        isleft,isbottom = issubplotedge(rows,cols,i+1)
        at = AnchoredText(str(int(Ur[i]))+'$\lambda$',
                  prop=dict(size=12), frameon=False,
                  loc=1,
                  )
        ax.add_artist(at)
        detect_area = Rectangle(
            (n.log10(0.1),Urspecz.min()),
            n.log10(0.25)-n.log10(0.1),
            Urspecz.max()-Urspecz.min(),fill=False,lw=2,ec='0.5')
        pl.xlim([-2,1])
        ax.add_patch(detect_area)
        if not isleft: pl.yticks([])
        if not isbottom: pl.xticks([])
    pl.colorbar(cax = cax,orientation='vertical')
    pl.subplots_adjust(wspace=0,hspace=0)
    pl.figtext(0.5,0.05,'log k')
    pl.figtext(0.05,0.5,'redshift',rotation='vertical')
    pl.draw()
    pl.savefig(file+'.png')
print "\a"*3
#for i in n.arange(nur-1,-1,-1):
#    if len(Urgrid[i])==0: del(Urgrid[i])
##do some nice plotting
#
#rows = ceil(n.sqrt(len(Urgrid)))
#cols = ceil(len(Urgrid)/rows)
#nplots = rows*cols
#pl.suptitle(', '.join(files)+'\n z = %4.1f'%z)
#for ri,ur in enumerate(Urgrid):
#    pl.subplot(rows,cols,ri+1)
#    for nps,Ps in enumerate(ur):
#        kparr = n.array([eta2kparr(d,z) for d in Urgrid_dly[ri][nps]])
#        pl.loglog(kparr,Ps,'k',alpha=0.2)
#    pl.loglog(kparr,n.mean(ur,axis=0),label=str(int(Ur[ri]))+'$\lambda$')
#    pl.legend()
#    pl.ylim([psmin,psmax])
#    isleft,isbottom = issubplotedge(rows,cols,ri+1)
#    if not isleft: pl.yticks([])
#    if not isbottom: pl.xticks([])
#pl.subplots_adjust(wspace=0,hspace=0)



#get u,v,ur coordinates for each 
    #sort power spectra into radial bins
    #plot each uv cell in ur subplots (we're looking for outliers)
