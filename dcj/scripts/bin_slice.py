#! /usr/bin/env python
import numpy as n, pylab as p,ephem
import sys, aipy as a,re
from pylab import *
import matplotlib.lines as lines
import itertools
from matplotlib.pyplot import cm
import matplotlib.animation as animation
r2h = 12/n.pi
i,j = 64,49
mybl=str(a.miriad.ij2bl(i,j))
jdjump = 2456942.
print mybl
def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])
def hms(h):
    return ephem.hours(h)
def dB(x):
    return 10*n.log10(x)
data = {}
time = {} #time is the lst axis... sorry
jd_times = {}
jds = []
for filename in sys.argv[1:]:
    jd = int(n.round(file2jd(filename)))
    t = file2jd(filename)
#    if jd<jdjump:continue  #focus on my epoch
#    print 'Reading', filename
    try:
        npz = n.load(filename)
    except:
        print "    failed to load"
    try:
        npz[mybl]

    except(KeyError):
        print "baseline {i}_{j} not found in {f}, skipping".format(i=i,j=j,f=filename)
        continue
    try:
        ntimes = len(npz[mybl])
        print npz['t'+mybl]
        dt = n.diff(npz['t'+mybl])[0]/n.pi #convert radians to fraction of a day
        data[jd].append(npz[mybl])
        time[jd].append(npz['t'+mybl])
        jd_times[jd].append(n.arange(0,ntimes)*dt + file2jd(filename))

    except(KeyError):
	try:
            data[jd] = [npz[mybl]]
            time[jd] = [npz['t'+mybl]]
            jd_times[jd] = [n.arange(0,ntimes)*dt + file2jd(filename)]
        except(KeyError):
            continue
if len(data)==0:
    print "ERROR: no data found"
    sys.exit()
for jd in data:
    data[jd] = n.concatenate(data[jd])#turna list of arrays into a single array
    time[jd] = n.concatenate(time[jd])
    jd_times[jd] = n.concatenate(jd_times[jd])
jd_samples = n.concatenate([jd_times[jd] for jd in jd_times])
#####
## GRID THE DATA
#make my lst grid    
lsts = n.concatenate([time[jd] for jd in time])
lstmin,lstmax,dlst = lsts.min(),lsts.max(),n.min(n.abs(n.diff(lsts)))
lst_grid = n.arange(lstmin*0.95,lstmax*1.05,dlst*2) #note I am tweaking the grid size a little bit after looking at the lst x day plot below

#make a jd grid
jds = n.array(time.keys())
jd_grid = n.arange(jds.min()-1,jds.max()+2,1)
daynums = jd_grid-jd_grid.min()
#grid the data into lst,jd bins
gridded_data = n.zeros((lst_grid.size,jd_grid.size)).astype(n.complex64) #dimensions lstbins x jds
counts = n.zeros_like(gridded_data).astype(int)
jd_inds = n.digitize(jds,jd_grid)
for i,jd in enumerate(jds):
    lst_inds = n.digitize(time[jd],lst_grid)
    gridded_data[lst_inds,jd_inds[i]] += data[jd]
    counts[lst_inds,jd_inds[i]] += 1
gridded_data[counts>0] /= counts[counts>0]
gridded_data = n.ma.masked_where(gridded_data==0,gridded_data)
jds = n.array(time.keys())

##### DATA IS GRIDDED


### MAKE A MEDIAN MODEL
lst_model = n.ma.median(gridded_data,axis=1)
lst_model.shape += (1,)
lst_model = n.repeat(lst_model,gridded_data.shape[1],axis=1)
res = gridded_data-lst_model
lst_model_errs = n.std(res,axis=1)
figure()
errorbar(lst_grid*r2h,lst_model[:,0],yerr=lst_model_errs)
xlabel('lst')
ylabel('Jy')
#### PLOT THE GRID
figure()
subplot(131)
imshow(n.abs(gridded_data),aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
colorbar()
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
#hist(n.abs(gridded_data.ravel()),bins=100)
subplot(132)
imshow(gridded_data.real,aspect='auto',interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
colorbar()
subplot(133)
imshow(n.abs(gridded_data-lst_model),aspect='auto',interpolation='nearest',
    extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h),vmax=0.05)
colorbar()


### Make a 2D histogram in lst
linear_data = n.concatenate([data[jd] for jd in jds])
linear_data = n.ma.masked_where(linear_data==0,linear_data)
Jy_range = n.array([-3,3])*n.std(linear_data.real)
Jy_res = n.diff(Jy_range)/25
lst_range = [lst_grid.min(),lst_grid.max()]
Jy_grid = n.arange(Jy_range[0],Jy_range[1],Jy_res)
H,lst_edges,Jy_edges = n.histogram2d(lsts[linear_data!=0],linear_data[linear_data!=0].real,
        bins=[lst_grid,Jy_grid])
##Make a basic plot of the slice vs lst and JD
figure()
subplot(221)
plot(lsts,linear_data,'.')
xlabel('lst')
subplot(222)
plot(jd_samples,n.concatenate([data[jd] for jd in jds]),'.')
xlabel('jd')
subplot(313)
imshow(H.T,aspect='auto',interpolation='nearest',extent=(lst_edges.min()*r2h,lst_edges.max()*r2h,Jy_edges.min(),Jy_edges.max()),cmap='cubehelix')
xlabel('lst [h]')
colorbar()
show()
sys.exit()

### Flag the residuals
ncut = 3
if False:
    ###Flag vs jd
    v = n.sqrt(n.mean(res*n.conj(res),axis=1))
    v.shape += (1,)
    v = n.repeat(v,res.shape[1],axis=1)
    res_m = n.ma.masked_where(n.logical_or(n.abs(res.real)>(ncut*v/n.sqrt(2)),
                            n.abs(res.imag)>(ncut*v/n.sqrt(2))),res)
if False:
    ### Flag vs lst
    v = n.sqrt(n.mean(res*n.conj(res),axis=0))
    v.shape = (1,)+v.shape
    v = n.repeat(v,res.shape[0],axis=0)
    res_m = n.ma.masked_where(n.logical_or(n.abs(res.real)>(ncut*v/n.sqrt(2)),
                            n.abs(res.imag)>(ncut*v/n.sqrt(2))),res)


### Flag using lst and jd
if True:
    v = n.sqrt(n.mean(res*n.conj(res)))
    res_m = n.ma.masked_where(n.abs(res)>(ncut*v),res)
print "flagging fraction: {frac}".format(frac=n.sum(res_m.mask-res.mask)/float(res_m.size))


### Plots of the residual.  There seems to be a slow phase drift.
if False:
    figure()
    subplot(131)
    imshow(res.real,aspect='auto',interpolation='nearest',vmin=-0.01,vmax=0.01)
    subplot(132)
    imshow(res_m.mask-res.mask,aspect='auto',cmap='gray_r',interpolation='nearest')
    subplot(133)
    imshow(res_m.real,aspect='auto',interpolation='nearest',vmin=-0.01,vmax=0.01)
    colorbar()
    
    figure()
    subplot(131)
    imshow(n.ma.abs(res_m),aspect='auto',interpolation='nearest',vmax=0.01)
    
    subplot(132)
    imshow(n.ma.angle(res_m),aspect='auto',interpolation='nearest')
    colorbar()
    subplot(133)
    res_m_frt = n.fft.fft(res_m,axis=0)
    imshow(n.flipud(n.abs(res_m_frt)),aspect='auto',interpolation='nearest',
    extent=(daynums.min(),daynums.max(),0,1/n.diff(lst_grid*r2h/60)[0]/2))
    print n.diff(lst_grid*r2h/60)[0],1/n.diff(lst_grid*r2h/60)[0]/2
    colorbar()

if False: 
    ### compare residuals at different times
    figure()
    #looks like 90 degrees from beginning to end
    plot(lst_grid*r2h,lst_model[:,0])
    plot(lst_grid*r2h,gridded_data[:,10],'g')
    plot(lst_grid*r2h,gridded_data[:,57],'r')
    plot(lst_grid*r2h,n.ma.mean(res_m[:,:10],axis=1),'--g')
    plot(lst_grid*r2h,n.ma.mean(res_m[:,56:],axis=1),'--r')
#its an ever so slight gain shift...
G = dB(n.ma.abs(gridded_data/lst_model))
G.mask[n.ma.abs(G)>5] = 1 #flag gains over 5dB and under -5dB
G_model_parms = n.ma.polyfit(daynums,n.mean(G,axis=0),1)
print "Gain model slope (dB/day):",G_model_parms[0]
G_model = n.poly1d(G_model_parms)
figure()
subplot(211)
imshow(G,aspect='auto',vmin=-1,vmax=1)
colorbar()
subplot(212)
gain_vs_day = n.ma.mean(G,axis=0)
plot(daynums,gain_vs_day)
plot(daynums,G_model(daynums),'k')
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
xlabel('daynum')

figure()
### plot the gain vs temps
try:
    temps = n.loadtxt('../Temperatures/nightly_temps_2456969_2456929.txt')
    plot(jd_grid,gain_vs_day,label='gain')
    #twinx()
    plot(temps[:,0],(temps[:,1]-n.mean(temps[:,1]))*(-0.06),'.k',label='nightly temp min')

except:
    print "couldn't find temps! skipping "


###gain corrected LST model
res_gain_corrected = gridded_data/10**(G_model(daynums)/10) - lst_model
res_gain_corrected.mask = res_m.mask
figure()
subplot(121)
imshow(res_m.real,aspect='auto',vmax=0.01)
subplot(122)
imshow(res_gain_corrected.real,aspect='auto',vmax=0.01)
colorbar()

### Check for xtalk
figure()
errorbar(daynums,n.mean(res_gain_corrected,axis=0).real,
            yerr=n.std(res_gain_corrected.real,axis=0),label='real')
errorbar(daynums,n.mean(res_gain_corrected,axis=0).imag,
            yerr=n.std(res_gain_corrected.imag,axis=0),label='imag')

if True:
    ### a basic test of integrating down
    figure()
    subplot(211)
    plot(lst_grid*r2h,n.ma.mean(n.abs(res.real),axis=1),label='inst')
    plot(lst_grid*r2h,n.abs(n.ma.mean(res,axis=1).real),label='all time')
    N_vs_lst = n.sum(n.logical_not(res.mask),axis=1)
    #plot(lst_grid*r2h,n.ones_like(lst_grid)*n.mean(n.abs(res),axis=1)/n.sqrt(len(jd_grid)),'--')
    plot(lst_grid*r2h,n.ma.mean(n.ma.abs(res.real),axis=1)/n.sqrt(N_vs_lst),'--k',label='2sigma')
    subplot(212)
    plot(lst_grid*r2h,n.ma.mean(n.abs(res_m.real),axis=1))
    plot(lst_grid*r2h,n.abs(n.ma.mean(res_m,axis=1).real))
    N_vs_lst = n.sum(n.logical_not(res_m.mask),axis=1)
    plot(lst_grid*r2h,n.mean(n.abs(res_m.real),axis=1)/n.sqrt(N_vs_lst),'--k',label='2sigma')
    xlabel('lst')
    show()
    sys.exit()
if False:
    ### Check that the residuals integrate down like noise
    nchoose = 10 #for each time length, draw and average this many samples of the residuals
    res = gridded_data - lst_model
    mylst = 4.4/r2h
    mylst_i = n.abs(lst_grid-mylst).argmin()
    integrated_residual = n.zeros((len(jds),len(lst_grid))).astype(n.complex64)
    res_lst_vs_time = [] #list of lsts each having rms vs number of days
    res_lst_vs_time_err =[]
    for i in range(len(lst_grid)):
        res_vs_time = []
        res_vs_time_err = []
        for j in range(len(jds)):
            vals = n.array([n.random.choice(res[i,:],j) for m in xrange(nchoose)])#sample,time
#            vals = res[i,:j]
            res_vs_time.append(n.mean(vals))
#            res_vs_time_err.append(n.std(n.mean(vals,axis=1)))
#            if i==mylst_i:
#                print res_vs_time[-1]
        res_lst_vs_time.append(res_vs_time)
        res_lst_vs_time_err.append(res_vs_time_err)
    res_lst_vs_time = n.array(res_lst_vs_time)
    figure()
    subplot(211)
    plot(res_lst_vs_time[mylst_i,:])
    xlabel('number of nights')
    ylabel('mean residual')
    subplot(212)
    loglog(n.abs(res_lst_vs_time.T))
    show()
    sys.exit()


#lets look at the day to day correlation
# I want to see a matrix of day_n x day_n, averaged over lst...
def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()
def cov2(m1,m2):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X1 = n.array(m1, ndmin=2, dtype=n.complex)
    X2 = n.array(m2, ndmin=2, dtype=n.complex)
    X1 -= X1.mean(axis=1)[(slice(None),n.newaxis)]
    X2 -= X2.mean(axis=1)[(slice(None),n.newaxis)]
    N = X1.shape[1]
    fact = float(N - 1)
    return (n.dot(X1, X2.T.conj()) / fact).squeeze()

C = cov(gridded_data.T)
print gridded_data.shape,C.shape
figure()
subplot(121)
imshow(n.abs(C),interpolation='nearest')
res_C = cov(gridded_data-lst_model)
subplot(122)
imshow(n.abs(res_C),interpolation='nearest')
figure()
d = (gridded_data-lst_model).ravel().real
d = d[d!=0]
H,X = histogram(d,bins=100)
H_max = H.max()
X = X[1:]
MEAN = n.sum(X*H)/n.sum(H)
WIDTH = n.sqrt(n.abs(n.sum((X-MEAN)**2*H)/n.sum(H)))
WHM = (X[H>H_max/2].max()-X[H>H_max/2].min())/2
MAX = n.mean(H[n.abs(X)<1e-3])
print H.max(),H[n.abs(X)<1e-3]
fit = lambda x: MAX*n.exp(-(x-MEAN)**2/(2*WIDTH**2))
fitnarrow = lambda x: MAX*n.exp(-(x-MEAN)**2/(2*WHM**2))
#hist(d[d!=0],bins=100,histtype='step')
semilogy(X,H)
plot(X,fit(X))
#plot(X,fitnarrow(X))
ylim([H[H>0].min(),H.max()])

#### PLOT AN ANIMATION
#make a Zheng-like thingy
fig,ax = subplots()
color=iter(cm.rainbow(np.linspace(0,1,num=len(jd_grid))))
line, = ax.plot(n.real(gridded_data[0,:]),n.imag(gridded_data[0,:]),'x')

if True:
    from scipy.interpolate import interp1d
    # smooth out each day
    smoothed_data = []
    smooth_factor = 10
    sm_lsts = n.linspace(lst_grid.min(),lst_grid.max(),num=len(lst_grid)*smooth_factor)
    for i in xrange(len(jd_grid)):
        d = gridded_data.data[:,i]
        if n.sum(d!=0)<2:
            smoothed_data.append(n.zeros_like(sm_lsts))
            continue
        #trim off edge zeros
        left,right = (d.size-len(n.trim_zeros(d,'f'))),(d.size-len(n.trim_zeros(d,'b')))
        d = d[left:-right]
        l = lst_grid[left:-right]
#        real_model = interp1d(lst_grid,gridded_data[:,i].real,kind='cubic')
#        imag_model = interp1d(lst_grid,gridded_data[:,i].imag,kind='cubic')
        real_model = interp1d(l,d.real,kind='cubic',bounds_error=False)
        imag_model = interp1d(l,d.imag,kind='cubic',bounds_error=False)
        smoothed_data.append(real_model(sm_lsts)+1j*imag_model(sm_lsts))
    smoothed_data = n.vstack(smoothed_data).T
    smoothed_data = n.ma.masked_invalid(smoothed_data)
    #smoothed_data = n.ma.masked_where(n.abs(smoothed_data)<1e-3,smoothed_data)
if False:
    figure()
    imshow(n.abs(smoothed_data),aspect='auto')
    colorbar()
    show()
    sys.exit()
if False:
    figure()
    plot(sm_lsts,smoothed_data[:,4].real,'.')
    plot(lst_grid,gridded_data[:,4].data.real)
    plot(sm_lsts,smoothed_data[:,4].imag,'.')
    plot(lst_grid,gridded_data[:,4].data.imag)
    print smoothed_data[:,4]
    show()
    sys.exit()
if True:
#show each days trajectory as a line of a different color
# use smoothed_data
    Nrecent = 10*smooth_factor
    d = smoothed_data[:Nrecent,:]
    linetrails = []
    for lst in xrange(d.shape[1]):
        c=next(color)
        linetrails += plot(d[:,i].real,d[:,i].imag,color=c)
    xlim([smoothed_data.real.min(),smoothed_data.real.max()])
    ylim([smoothed_data.imag.min(),smoothed_data.imag.max()])
    print smoothed_data.real.max(),smoothed_data.imag.max(),n.abs(smoothed_data).max()
#    scatter(smoothed_data.real.ravel(),smoothed_data.imag.ravel(),marker='.')
    def trails_animate(lst_i):
        d=smoothed_data[:lst_i,:]
        if d.shape[0]>Nrecent:
            d = d[-1*Nrecent:,:]
        for i in xrange(d.shape[1]): 
            linetrails[i].set_xdata(d[:,i].real)
            linetrails[i].set_ydata(d[:,i].imag)
        ax.set_title(hms(sm_lsts[lst_i]))


if False:
#show each days trajectory as a line of a different color
    Nrecent = 10
    d = gridded_data[:Nrecent,:]
    linetrails = []
    for lst in xrange(d.shape[1]):
        c=next(color)
        linetrails += plot(d[:,i].real,d[:,i].imag,color=c)
    def trails_animate(lst_i):
        d=gridded_data[:lst_i,:]
        if d.shape[0]>Nrecent:
            d = d[-1*Nrecent:,:]
        for i in xrange(d.shape[1]): 
            linetrails[i].set_xdata(d[:,i].real)
            linetrails[i].set_ydata(d[:,i].imag)
        

if False:
#this option prints with decaying alpha
    Nrecent = 10  
    d=gridded_data[:Nrecent,:]
    print d.shape
    trails = []
    for i in xrange(d.shape[0]):
        trails+=plot(d[i,:].real,d[i,:].imag,'.',c='b')#,alpha=n.linspace(0,1,num=d.shape[0]))
    print trails
    print Nrecent,len(trails)   
    alphas = n.linspace(0.5,1,num=Nrecent)
    #todo just plot each step seperately, delete when older than Nrecent, alpha decays with time.
    def trails_animate(lst_i):
        #get the most recent points
        #plot as dots them with alpha increasing into the past
        d=gridded_data[:lst_i,:]
        if d.shape[0]>Nrecent:
            d = d[-1*Nrecent:,:]
        for i in xrange(d.shape[0]): #for each lst, plot the jds with an increasing alpha
            trails[i].set_xdata(d[i,:].real)
            trails[i].set_ydata(d[i,:].imag)
            trails[i].set_alpha(alphas[i])
        return trails,
        
def animate(lst_i):
    line.set_xdata(n.real(gridded_data[lst_i,:]))
    line.set_ydata(n.imag(gridded_data[lst_i,:]))
    ax.set_title(hms(lst_grid[lst_i]))
    return line,
def init():
    line.set_ydata(n.real(n.ma.array(gridded_data[0,:],mask=True)))
    line.set_xdata(n.imag(n.ma.array(gridded_data[0,:],mask=True)))
    return line,
#ani = animation.FuncAnimation(fig,trails_animate,xrange(Nrecent,len(lst_grid)),interval=10,blit=False,init_func=init,repeat=True)
ani = animation.FuncAnimation(fig,trails_animate,xrange(Nrecent,len(lst_grid)*smooth_factor),interval=10,blit=False,init_func=init,repeat=True)

if False:
    animation_file = 'vis_grid_anim_{jdmin}_to_{jdmax}_lines.gif'.format(jdmin=int(jd_grid.min()),
                                                        jdmax=int(jd_grid.max()))
    print "saving animation to:",animation_file
    ani.save(animation_file,writer='imagemagick',fps=10)
#for i,lst in enumerate(lst_grid):
#    markers = itertools.cycle(lines.Line2D.markers.keys())
#    marker = markers.next()
#    c = next(color)
#    scatter(n.real(gridded_data[i,:]),n.imag(gridded_data[i,:]),c=c)#,marker=marker)
#    hist(n.real(gridded_data[i,:]),bins=10,histtype='step',color=c)
print "finished"
show()
