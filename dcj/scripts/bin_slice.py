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

data = {}
time = {}
jds = []
for filename in sys.argv[1:]:
    jd = int(n.round(file2jd(filename)))
    if jd<jdjump:continue  #focus on my epoch
#    print 'Reading', filename
    try:
        npz = n.load(filename)
    except:
        print "    failed to load"
    try:
        data[jd].append(npz[mybl])
        time[jd].append(npz['t'+mybl])
    except(KeyError):
	try:
            data[jd] = [npz[mybl]]
            time[jd] = [npz['t'+mybl]]
        except(KeyError):
            continue
for jd in data:
    data[jd] = n.concatenate(data[jd])#turna list of arrays into a single array
    time[jd] = n.concatenate(time[jd])


#####
## GRID THE DATA
#make my lst grid    
lsts = n.concatenate([time[jd] for jd in time])
lstmin,lstmax,dlst = lsts.min(),lsts.max(),n.min(n.abs(n.diff(lsts)))
lst_grid = n.arange(lstmin*0.95,lstmax*1.05,dlst*2) #note I am tweaking the grid size a little bit after looking at the lst x day plot below

#make a jd grid
jds = n.array(time.keys())
jd_grid = n.arange(jds.min()-1,jds.max()+2,1)
print jd_grid.min(),jds.min(),jds.max(),jd_grid.max()
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

#### PLOT THE GRID
subplot(131)
imshow(n.abs(gridded_data),aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
#hist(n.abs(gridded_data.ravel()),bins=100)
ax2 = subplot(132)
imshow(gridded_data.real,aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
ax2.set_yticklabels([])
subplot(133)
imshow(n.abs(gridded_data-lst_model),aspect='auto',interpolation='nearest',
    extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h),vmax=0.05)
colorbar()
subplots_adjust(hspace=0)

### PLOT a lst slice.
mylst = 4.36/r2h #hours of lst
mylst_i = n.abs(lst_grid-mylst).argmin()
print mylst_i
figure()
for i in range(mylst_i-3,mylst_i+3):
    #lstname = hms(lst_grid[i])
    lstname = str(n.round(lst_grid[i]*r2h,2))
    plot(jd_grid-jd_grid.min(),n.abs(gridded_data-lst_model)[i,:],label=lstname)
legend()
text(0.92,-0.07,"+%i"%jd_grid.min(),fontsize=10,transform=gca().transAxes)
title('lst slices')
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
def corrcoef(x):
    C = cov(x)
    return C/n.sqrt(n.outer(n.diag(C),n.diag(C)))
### PLot COVARIANCES
C = cov(gridded_data.T)
print gridded_data.shape,C.shape
figure()
subplot(131)
imshow(n.abs(C),interpolation='nearest')
res_C = cov(gridded_data-lst_model)
subplot(132)
imshow(n.abs(res_C),interpolation='nearest')
subplot(133)
corr = corrcoef(gridded_data.T)
imshow(n.abs(corr),interpolation='nearest')

### plot/fit a historgram to the residuals
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
plot(X,H,label='res')
plot(X,fit(X),label='fit')
legend()
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
