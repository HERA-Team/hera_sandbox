#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import aipy as a, numpy as n, pylab as p, sys, optparse, glob, ipdb, ephem, capo as C
o=optparse.OptionParser()
o.set_usage("plot_tsys.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, cmap=True,chan=True)
o.add_option('--plot',dest='plot',default=False, action='store_true',\
    help='Outputs plot to X-Window and saves plot')
o.add_option('--output', type='string', default='',
    help='output directory for image files (default "")')
o.add_option('--freqs',dest='freqs',default=False, action='store_true',\
    help='Adds Frequency Axis to Tsys Plot')
o.add_option('--vline',default=False,action='store_true',\
    help='Emphasizes chosen channel range')
opts,args=o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal,.1,.1,1)



FREQ_AXIS=opts.freqs
VLINE=opts.vline


freqs = None
#dsets=glob.glob('/data3/PAPER/psa64/lstbin_omnical_2/lstbinX0/sep0,1/*.uvAL')
#dsets = n.sort(dsets)
dsets=args
times=[]

plot_n={'jd':[], 'lst':[], 'lst_cal':[], 'cnt':[], 'nchan':[],'var':[]}
plot_x={}

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
            plot_n['var'].append(uv['var'])

lsts,cnt,var = n.array(plot_n['lst']), n.array(plot_n['cnt']), n.array(plot_n['var'])
chans=n.arange(plot_n['nchan'][0])
lsts = n.where(lsts > 5, lsts - 2*n.pi, lsts) * 12 / n.pi
fqs = n.linspace(.1,.2,var.shape[1])
jy2T = C.pspec.jy2T(fqs)
jy2T.shape = (1,203)
extent = (1e3*fqs[0], 1e3*fqs[-1], lsts[-1], lsts[0])
Tlst = n.sqrt(var)*jy2T
#three_sig_correction = 0.974
three_sig_correction = 1/1.34
Tsys = Tlst * 1e-3 * n.sqrt(43 * 100. / 203 * 1e6) / three_sig_correction
Trms = Tlst / n.sqrt(cnt.clip(1,n.Inf))

p.figure(figsize=(7,4.2))
chunk = 83
for i in xrange(40,600,chunk):
    print i, n.average(lsts[i:i+chunk])
    lst_hr = int(n.around(n.average(lsts[i:i+chunk])))
    Tsys_avg = n.sqrt(n.average(Tsys[i:i+chunk]**2, axis=0))
    valid = n.where(Tsys_avg > 213, 1, 0)
    #p.plot(1e3*fqs.compress(valid), Tsys_avg.compress(valid), label='LST = %02dh00'%lst_hr)
    p.plot(chans.compress(valid), Tsys_avg.compress(valid), label='LST = %02dh00'%lst_hr)

if VLINE:
    nums=opts.chan.split('_')
    vlines=[float(v) for v in nums]
    p.vlines(vlines,200,2000,linestyles='dashed')
p.legend(loc='best', fontsize='medium')
p.xlabel(r'${\rm Channel}$', fontsize=16)
if FREQ_AXIS:
    ax=p.twiny()
    ax.plot(freqs,n.ones(len(freqs*1e3)))
    ax.cla()
    ax.set_xlabel(r'${\rm Frequency\ [MHz]}$', fontsize=16)
#p.xlim(125,175)
p.ylim(200,2000)
p.grid()
p.ylabel(r'${\rm T}_{\rm sys}$', fontsize=16)
p.subplots_adjust(0.12, 0.15, .95, .85)

fig_file='Tsys_chan_'+opts.chan
if not opts.output == '':
    fig_file= opts.output + '/' + fig_file
p.savefig(fig_file)
#p.show()

Tsys_jy = Tsys / 1e-3 / jy2T

#n.savez('tsys_model.npz', tsys=Tsys, freq=fqs, lsts=lsts)
#n.savez('tsys_model_jy.npz', tsys_jy=Tsys_jy, freq=fqs, lsts=lsts)
if False:
    p.subplot(121)
    C.arp.waterfall(Tsys, mx=700, drng=500, mode='lin', extent=extent)
    p.xlim(110,185)
    p.xlabel('Frequency [MHz]')
    p.xticks([125, 150, 175])
    p.colorbar()

    p.subplot(122)
    C.arp.waterfall(Trms, mx=30, mode='lin', extent=extent)
    p.xlim(110,185)
    p.xlabel('Frequency [MHz]')
    p.xticks([125, 150, 175])
    p.colorbar()
    p.show()
