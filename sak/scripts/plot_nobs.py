#! /usr/bin/env python
"""

Creates waterfall plot of N_obs from Miriad UV files.


"""


import aipy as a, numpy as n, sys, optparse, glob, ephem
from matplotlib import pylab as p
o=optparse.OptionParser()
o.set_usage("plot_nobs.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, cmap=True)
o.add_option('--plot',dest='plot',default=False, action='store_true',\
    help='Outputs plot to X-Window and saves plot')
o.add_option('--var',dest='var',default=False, action='store_true',\
    help='Plot in variance instead of occupancy')
opts,args=o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal,.1,.1,1)
freqs = None
#dsets=glob.glob('/data3/PAPER/psa64/lstbin_omnical_2/lstbinX0/sep0,1/*.uvAL')
#dsets = n.sort(dsets)
dsets=args
times=[]

plot_n={'jd':[], 'lst':[], 'lst_cal':[], 'cnt':[], 'nchan':[]}
plot_x={}

for F in dsets:
    print 'Reading', F.split('/')[-1]
    uv = a.miriad.UV(F)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    if freqs is None:
        aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
        freqs = aa.get_afreqs()
        freqs *= 1e3 #convert to MHz
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
            if opts.var: plot_n['cnt'].append(n.log10(uv['var']))
            else: plot_n['cnt'].append(uv['cnt'])
            plot_n['nchan'].append(uv['nchan'])

if len(times) ==0:
        print('No data to plot.')
        sys.exit(0)
nchan = plot_n['nchan'][0]

cmap= p.get_cmap(opts.cmap)
fig, ax1= p.subplots()
chans=n.arange( max( n.array( plot_n['nchan']  ) ) )  
n_color=n.array(plot_n['cnt']).max()
levels=n.arange(n_color)
cnt_plot =  n.array(plot_n['cnt'])
#cnt_plot=n.reshape(plot_n['cnt'],(len(chans),len(plot_n['jd'])))
#img=ax1.contour( chans , plot_n['jd'], cnt_plot.T , levels )
lsts = n.array(plot_n['lst'])
img=ax1.imshow(cnt_plot,aspect='auto',interpolation='nearest',cmap=cmap,
    extent=(0,203,lsts.min()*12/n.pi,lsts.max()*12/n.pi))
#ax1.set_yticklabels(str(plot_n['lst']))
#str_ticks=[ str(plot_n['lst'][k]) for k in xrange(len(times))]
#ax1.set_yticklabels(str_ticks)
p.gca().invert_yaxis()
if opts.var: ax1.set_title('Observation log10(Variance)')
else: ax1.set_title('Observation Counts')
ax1.set_ylabel('LST')
ax1.set_xlabel('Frequency [MHz]')
fig.subplots_adjust(right=0.85,left=.15)
cbar_ax=fig.add_axes([0.9,0.25,.02,.5])
fig.colorbar(img,cax=cbar_ax)
if opts.var: fig.savefig('N_var_waterfall.png',format='png')
else: fig.savefig('N_obs_waterfall.png',format='png')

p.figure()
p.plot(lsts*12/n.pi,cnt_plot[:,100],'.')
p.xlabel('lst')
p.ylabel('count')
#p.gca().invert_yaxis()

# make a pcolor plot
if False: #experimenting with trying to map lst to pixels nicely
    FREQS,LSTS = n.meshgrid(freqs,lsts*12/n.pi)
    p.figure()
    p.pcolor(FREQS,LSTS,cnt_plot)

if opts.plot:
    p.show()
p.close()
