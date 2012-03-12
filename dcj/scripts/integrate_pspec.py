#test of integration of power spectrum
#inputs: list of casa pspec files, list of kparr ranges with radial bin setup. (eg [0.1,0.24,25,9.5] would be
#kparrlow,kparrhigh,ur,z)
#eg draw boxes in the z x k plot and integrate down. 

#load a file
#peform the usual k computation and radial binning
#loop over the bins I am interested in and grab those samples
#for each bin, get a list of samples and times
#save and plot somehow

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
intlog = {}
intbins = [[0.1,0.25,40,9],[0.1,0.25,20,9],[0.1,0.25,40,8],[0.4,0.7,40,9],[1,1.5,20,9],[1,1.5,40,8],[10**0.4,10**0.8,30,8]]

for file in files:
    #load file
    print file
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

    Zs = n.array([f212z(f) for f in fs])
    for band,f in enumerate(fs):
    #    f = csys.toworld([0,0,band,0])['numeric'][2]
        z = f212z(f)
        binzis = [plop(b[3],Zs) for b in intbins]
#        if n.sum(n.abs(z-binzs)>0.25)>0: continue
        #find the delay axis
        dly = n.array([csys.toworld([0,0,band,i])['numeric'][3] for i in\
            range(uveta.shape[3])])
        #build kparr,z coordinate arrays for display purposes
        kparr = n.array([eta2kparr(d,z) for d in dly])
        

        for i,UK in enumerate(UrspecK):
            k = n.sqrt(kparr**2 + u2kperp(Ur[i],z)**2) 
            if dologk:
                logk  = n.logspace(n.log10(k[k>0].min()),
                            n.log10(k.max()),
                            num=nkparr)
                UK[band,:] = logk
            else:
                UK[band,:] = k
        Urspecz[:,band,:] = z

#        if not (band in binzis): continue
        print f,z;flush()
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
            for bini,bin in enumerate(intbins):
                k_min,k_max,ur,intz = bin
                kp_min = n.sqrt(k_min**2 - u2kperp(ur,intz)**2)
                kp_min = n.sqrt(k_max**2 - u2kperp(ur,intz)**2)
                ui = plop(Ur,ur)
                intband = plop(Zs,intz)
                k = UrspecK[ui,band,:]
#                Urspec_snapshot = Urspec.copy()
#                Urspec_snapshot[Urspecn>0] /= Urspecn[Urspecn>0]
#                PS = Urspec_snapshot[plop(Ur,length(UVworld)),intband,:] 
               # del(Urspec_snapshot)
                if plop(Ur,length(UVworld))==ui and band==intband:
                    kpi_min = plop(k,k_min)
                    kpi_max = plop(k,k_max)
                    N = Urspecn[plop(Ur,length(UVworld)),intband,0]
                    intlog[bini] = intlog.get(bini,[]) + [[t,n.mean(Urspec[plop(Ur,length(UVworld)),intband,kpi_min:kpi_max+1])/N,N,n.mean(ps[kpi_min:kpi_max+1])]]

            
        if uveta.max()>psmax:psmax=uveta.max()
        if uveta.min()<psmin:psmin=uveta.min()
Urspec[Urspecn>0] /= Urspecn[Urspecn>0]


#    Urspec = Urspec.filled(0)
pl.figure(16)
pl.clf()
ax = pl.subplot(111)
bincolors = []
for bini,bin in enumerate(intbins):
    if not intlog.has_key(bini):continue
    bin = n.round(n.array(bin),2)
    S = n.array(intlog[bini])
#    pl.subplot(211)
    ax.loglog(n.arange(1,S.shape[0]+1),S[:,1],'.',label='_'.join(map(str,bin)))
    lastcolor=ax.lines[-1].get_color()
    bincolors.append(lastcolor)
#    pl.subplot(212)
#    pl.loglog(n.arange(1,S.shape[0]+1),S[:,1],'.',label='_'.join(map(str,bin)))
#    pl.loglog(n.arange(1,S.shape[0]+1),n.cumsum(S[:,3])/n.cumsum(n.arange(1,S.shape[0]+1)))
    pl.loglog(n.arange(1,S.shape[0]+1),n.cumsum(S[:,3])/n.arange(1,S.shape[0]+1),'-'+lastcolor)
#    pl.loglog(n.arange(1,S.shape[0]+1),S[0,1]/n.arange(1,S.shape[0]+1),'--'+lastcolor)
#    pl.loglog(n.arange(1,S.shape[0]+1),S[:,2],label='_'.join(map(str,bin)))
pl.grid()
#pl.subplot(211)
pl.legend()
print "\a"
rects = {}
dz = n.diff(Zs)[0]
if True:
#    print type(Urspec)
    for bini,bin in enumerate(intbins):
        k_min,k_max,ur,intz = bin 
        rects[bini] = Rectangle(
        (n.log10(k_min),intz-dz/2),n.log10(k_max)-n.log10(k_min),dz,
        fill=False,lw=2,ec=bincolors[bini]
        )
    pl.figure(42)
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
        for bini,bin in enumerate(intbins):
            k_min,k_max,ur,intz = bin 
            if plop(Ur,ur)==i:
                ax.add_patch(rects[bini])
        if not isleft: pl.yticks([])
        if not isbottom: pl.xticks([])
    pl.colorbar(cax = cax,orientation='vertical')
    pl.subplots_adjust(wspace=0,hspace=0)
    pl.figtext(0.5,0.05,'log k')
    pl.figtext(0.05,0.5,'redshift',rotation='vertical')
    pl.draw()
    pl.savefig(file+'.png')
print "\a"*3
    
