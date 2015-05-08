#! /usr/bin/env python
"""

Creates waterfall plot of N_obs for LSTs with below average observations from Miriad UV files.


"""


import aipy as a, numpy as n, pylab as p, sys, optparse, glob, ipdb, ephem
from matplotlib.patches import Rectangle
o=optparse.OptionParser()
o.set_usage("plot_nobs.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, cmap=True)
o.add_option('--plot',dest='plot',default=False, action='store_true',\
    help='Outputs plot to X-Window and saves plot')
opts,args=o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal,.1,.1,1)
freqs = None
#dsets=glob.glob('/data3/PAPER/psa64/lstbin_omnical_2/lstbinX0/sep0,1/*.uvAL')
#dsets = n.sort(dsets)
dsets=args
times=[]

plot_n={'jd':[], 'lst':[], 'lst_cal':[], 'cnt':[], 'nchan':[]}
plot_x={}
flags={'zero':[],'nobs':[]}

for F in dsets:
    print 'Reading', F.split('/')[-1]
    uv = a.miriad.UV(F)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    if freqs is None:
        aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
        freqs = aa.get_afreqs()
        freqs *= 1e3 #convert to MHz
    #ipdb.set_trace()
    #
    for (uvw,t,(i,j)),d in uv.all():
        bl= '%d,%d,%d' % (i,j,uv['pol'])
        if len(times) ==0 or times[-1] != t:
            times.append(t)
            lst=uv['lst']
            aa.set_jultime(t)
            lst_cal=aa.sidereal_time()
            plot_n['jd'].append(t)
            plot_n['lst'].append(ephem.hours(uv['lst']))
            plot_n['lst_cal'].append(lst_cal)
            plot_n['cnt'].append(uv['cnt'])
            plot_n['nchan'].append(uv['nchan'])
            flags['zero'].append(uv['cnt'].tolist().count(0))
            flags['nobs'].append(n.sum(uv['cnt']))

if len(times) ==0:
        print('No data to plot.')
        sys.exit(0)
nchan = plot_n['nchan'][0]

lsts=n.array(plot_n['lst'])
print n.min(n.diff(lsts)),n.median(n.diff(lsts)),
dlst=n.median(n.diff(lsts))
print "total hours:",len(lsts)*dlst*12/n.pi
lst_grid=n.arange(0, 2*n.pi +dlst,dlst)
#lst_grid

lst_bins=n.digitize(lsts,lst_grid)

cnt_plot=n.zeros((len(lst_grid),nchan))
counts=n.zeros((len(lst_grid),nchan))
for j in xrange(nchan):
     cnt_plot[lst_bins,j] += n.array(plot_n['cnt'])[:,j]
     counts[lst_bins,j]+=1.
cnt_plot[counts>0] /= counts[counts>0]
print(counts.max())
#p.imshow(counts,aspect='auto')
#p.show()

cmap= p.get_cmap(opts.cmap)
fig, ax1= p.subplots()
chans=n.arange( max( n.array( plot_n['nchan']  ) ) )  
n_color=n.array(plot_n['cnt']).max()
levels=n.arange(n_color)
#cnt_plot=n.reshape(bad_zero_lst, (len(chans),num_zero) )
#cnt_plot=n.array(plot_n['cnt'])
img=ax1.imshow(cnt_plot,aspect='auto',interpolation='nearest',cmap=cmap,
    extent=(freqs.min(),freqs.max(),24,0),
    vmin=0.0,vmax=levels.max())
#ax1.set_yticklabels(str(plot_n['lst']))
#str_ticks=[ str(plot_n['lst'][k]) for k in xrange(len(times))]
#ax1.set_yticklabels(str_ticks)
ax1.set_title('Observation Counts')
ax1.set_ylabel('LST')
ax1.set_xlabel('Frequency [MHz]')
fig.subplots_adjust(right=0.85,left=.15)
cbar_ax=fig.add_axes([0.9,0.25,.02,.5])
fig.colorbar(img,cax=cbar_ax)

f_min=140
f_max=181

good_freqs=n.where(n.logical_and(freqs > f_min, freqs < f_max ) )[0]
thresholds=n.arange(0,levels.max()) +1
n_good=[]
for hold in thresholds:
    n_bad= n.array([ n.sum(cnt_plot[i,good_freqs] < hold ) for i in xrange(len(cnt_plot))])
    n_good.append(n.sum(n_bad ==0))

n_good.insert(0,len(lst_bins))
thresholds=n.concatenate([[0],thresholds])

n_good=n.array(n_good)*dlst*12./n.pi
p.figure()
p.plot(thresholds,n_good,'.')
p.xlabel('Threshold')
p.ylabel('Hours of Good LSTs')


vert_x=freqs[good_freqs[0]]
x_len=freqs[good_freqs[-1]] - freqs[good_freqs[0]]

vert_y=lst_grid.min()*24./(2*n.pi)
y_len=4*24./(2*n.pi)
ax1.add_patch(Rectangle( (vert_x,vert_y),x_len, y_len,alpha=1,facecolor='none',
            linewidth=1.5))

fig.savefig('N_obs_waterfall.png',format='png')


p.figure()
p.plot(n_bad,lst_grid*24./(2*n.pi),'.')
p.ylim([lst_grid.max()*24./(2*n.pi),lst_grid.min()*24./(2*n.pi)])
p.xlabel('N bad')
p.ylabel('lst')

ind_min= n.where( flags['zero'] == n.amax(flags['zero']))[0]
if True:
    p.figure()
    mean_cnt= n.mean(cnt_plot,axis=0)
    p.plot(freqs,mean_cnt/mean_cnt.max(),'-')
    for i in lst_bins[ind_min]:
        p.fill_between(freqs,cnt_plot[i,:]/n.max(cnt_plot[i,:]),alpha=.5)
    p.ylim([0,1.1])
    p.xlabel('Frequncy [MHz]')
    p.ylabel('count')
else:
    p.figure()
    p.plot(lsts*12/n.pi,n.mean(cnt_plot,axis=1),'.')
    p.ylim([-5,levels.max()+1])
    p.xlabel('lst')
    p.ylabel('count')

#ipdb.set_trace()
# make a pcolor plot
if False: #experimenting with trying to map lst to pixels nicely
    FREQS,LSTS = n.meshgrid(freqs,lsts*12/n.pi)
    p.figure()
    p.pcolor(FREQS,LSTS,cnt_plot)

if opts.plot:
    p.show()
p.close()
