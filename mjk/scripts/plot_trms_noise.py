#! /usr/bin/env python
"""

Creates waterfall plot of T_RMS from Simulated White Noise.


"""


import aipy as a, numpy as n, pylab as p, sys, optparse, glob, ipdb, ephem, capo
o=optparse.OptionParser()
o.set_usage("plot_trms.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, cmap=True,chan=True)
o.add_option('--plot',dest='plot',default=False, action='store_true',\
    help='Outputs plot to X-Window and saves plot')
o.add_option('--output', type='string', default='',
    help='output directory for image files (default "")')
o.add_option('--vline',default=False,action='store_true',\
    help='Emphasizes chosen channel range')
o.add_option('--band', default='0_202', action='store',
    help='Channels from which to Calculate full band Covariance')
o.add_option('--noise',  type="float", default=0,
    help='Creates White Noise spectrum. Input is rms in mk')
opts,args=o.parse_args(sys.argv[1:])


POL = 'I'
FRF_WIDTH=401
NOISE= opts.noise

try:
    rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []

def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, cnt, var, dat, flg = [], {}, {}, {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if bl in rmbls: continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: 
                lsts.append(lst)
                #var.append(uv['var'])
            if not dat.has_key(bl):
                 dat[bl],flg[bl],cnt[bl],var[bl] = [],[],[],[]
            dat[bl].append(d)
            flg[bl].append(f)
            cnt[bl].append(uv['cnt'])
            var[bl].append(uv['var'])
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    return n.array(lsts),cnt,var, dat, flg

def noise(size):
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

aa = a.cal.get_aa(opts.cal,.1,.1,1)
#dsets=glob.glob('/data3/PAPER/psa64/lstbin_omnical_2/lstbinX0/sep0,1/*.uvAL')
#dsets = n.sort(dsets)


freqs = None

dsets = {
    'even': [x for x in args if 'even' in x],
    'odd' : [x for x in args if 'odd' in x]
}

print 'Number of even data sets: {0:d}'.format(len(dsets['even']))
print 'Number of odd data sets: {0:d}'.format(len(dsets['odd']))
for dset_count in xrange(len(dsets['even'])):
        print dsets['even'][dset_count].split('/')[-1], dsets['odd'][dset_count].split('/')[-1]

uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
band_chans = a.scripting.parse_chans(opts.band, uv['nchan'])
inttime = uv['inttime'] * 4 # XXX hack for *E files that have inttime set incorrectly
#inttime= 8*60.
print 'inttime', inttime
del(uv)

afreqs = freqs.take(chans)
allfreqs = freqs.take(band_chans)
nchan_band = len(band_chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, allfreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(allfreqs)

antstr = 'cross'
lsts, cnt1, var1, data, flgs = {}, {},{},{},{}
days = dsets.keys()
for k in days:
    lsts[k],cnt1[k],var1[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)
    print data[k].keys()

if True:
#if False:
    # Align data sets in LST
    print [lsts[k][0] for k in days]
    lstmax = max([lsts[k][0] for k in days])
    for k in days:
        print k
        for i in xrange(len(lsts[k])):
            # allow for small numerical differences (which shouldn't exist!)
            if lsts[k][i] >= lstmax - .001: break
        lsts[k]= lsts[k][i:]
#cnt1[k], var1[k] , cnt1[k][i:], var1[k][i:]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = data[k][bl][i:],flgs[k][bl][i:]
            cnt1[k][bl], var1[k][bl] = cnt1[k][bl][i:],var1[k][bl][i:]
    print [len(lsts[k]) for k in days]
    j = min([len(lsts[k]) for k in days])
    for k in days:
        lsts[k]= lsts[k][:j]
#cnt1[k], var1[k] = , cnt1[k][:j], var1[k][:j]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = n.array(data[k][bl][:j]),n.array(flgs[k][bl][:j])
            cnt1[k][bl],var1[k][bl] = cnt1[k][bl][:j],var1[k][bl][:j]

lsts=lsts.values()[0]

if len(lsts) ==0:
        print('No data to plot.')
        sys.exit(0)

x, cnt, var = {},{},{}
if set(['even','odd']) == set(days):
    for bl in data['even']:
        d = n.copy(data['even'][bl][:,band_chans]-data['odd'][bl][:,band_chans])* jy2T
        c = n.copy(n.array(cnt1['even'][bl])[:,band_chans]) 
        v = n.copy(n.array(var1['even'][bl])[:,band_chans])
        if conj[bl]: d=n.conj(d)        
        x[bl] = n.transpose(d,[1,0])
        cnt[bl] = n.transpose(c,[1,0])
        var[bl] = n.transpose(v, [1,0])*(inttime/(8*60.)) ## Fix inttime for new FRF
bls_master=x.keys()
nbls = len(bls_master)
bls_master.sort()

wnx={}
if NOISE > 0:
        print 'Creating white noise at {0} mk rms'.format(NOISE)
        tmp_nx = {}
        bl1 = a.miriad.bl2ij(bls_master[0])
        beam_w_fr = capo.frf_conv.get_beam_w_fr(aa,bl1,ref_chan=0)
        t, firs, frbins, frspace = capo.frf_conv.get_fringe_rate_kernels(beam_w_fr, inttime, FRF_WIDTH)
        for k in days:
            wnx[k]={}
            tmp_nx[k] = {}
            for bl in bls_master:
                noise1 = noise(n.transpose(data[days[0]][bls_master[0]], [1,0]).shape)*6.6408 * NOISE ### mk 
                for count, ch in enumerate(band_chans):
                    noise1[count] = n.convolve(noise1[count], firs[count], mode='same')

                tmp_nx[k][bl] = noise1.copy() 
                wnx[k][bl]=noise1.copy()
    
        sub_noise={}
        for bl in bls_master:
                sub_noise[bl] = n.copy( wnx['even'][bl] - wnx['odd'][bl])
   
trms_data, trms_the = {}, {}
theo_temp= {}
noise_rms = {}
wnoise_rms = {}

nlst=len(lsts)
dlst= n.ceil(nlst/13.)

for bl in bls_master:
    trms_data[bl] = n.sqrt( n.mean(x[bl][band_chans,::dlst].conj() * x[bl][band_chans,::dlst],axis=1)/2. ).real
    noise_rms[bl] = n.sqrt( n.mean(sub_noise[bl][band_chans,::dlst].conj() * sub_noise[bl][band_chans,::dlst] ,axis=1)/2. ).real
 #   wnoise_rms[bl] =n.sqrt( n.mean(wnx['even'][bl][band_chans,::dlst].conj()*wnx['odd'][bl][band_chans,::dlst]  ,axis=1)).real
    trms_the[bl] = n.sqrt( n.mean(var[bl][band_chans,::dlst]/cnt[bl][band_chans,::dlst].clip(1,n.Inf)**2,axis=1) ) * jy2T[band_chans]
##GSM emission + antenna temp / sqrt(df*dt*cnt)
    theo_temp[bl]= (1.2e5*(allfreqs/.15)**(-2.8)+400)/(n.sqrt( 8*60. *100/203. *1e6))* n.mean(1./n.sqrt(cnt[bl][band_chans,::dlst].clip(1,n.Inf)),axis=1)


var_blavg = n.mean( [var[bl] for bl in bls_master],axis=0)
cnt_blavg = n.mean( [cnt[bl] for bl in bls_master],axis=0)

trms_the_blavg=n.sqrt( n.mean(var_blavg[band_chans,::dlst]/(cnt_blavg[band_chans,::dlst].clip(1,n.Inf)**2*nbls),axis=1)) * jy2T[band_chans]

diff_blavg = n.mean( [x[bl] for bl in bls_master],axis=0)
trms_blavg= n.sqrt( n.mean(diff_blavg.conj()[band_chans,::dlst]*diff_blavg[band_chans,::dlst] ,axis=1)/2.)

n_blavg = n.mean( [sub_noise[bl] for bl in bls_master],axis=0)
#e_blavg = n.mean( [wnx['even'][bl] for bl in bls_master],axis=0)
#o_blavg = n.mean( [wnx['odd'][bl] for bl in bls_master],axis=0)
noise_blavg= n.sqrt( n.mean(n_blavg.conj()[band_chans,::dlst]*n_blavg[band_chans,::dlst],axis=1)/2.)
#wnoise_blavg = n.sqrt(n.mean( e_blavg[band_chans,::dlst].conj()*o_blavg[band_chans, ::dlst],axis=1))
 
avg_trms= n.mean([ trms_data[bl] for bl in bls_master], axis=0)
gsm_data={}

for bl in bls_master:
    gsm_data[bl] = n.ma.masked_invalid(theo_temp[bl]/avg_trms)
    print 'baseline:', bl, 'mean, median, std ratio', n.mean(gsm_data[bl]), n.median(gsm_data[bl]), n.std(gsm_data[bl])

for bl in bls_master:
    print 'Ploting Trms for %d_%d'%a.miriad.bl2ij(bl)
    p.plot(band_chans, trms_data[bl], label='$T_{e-o}$')
    #p.plot(band_chans, noise_rms[bl],'-.', label='$T_{noise}$')
    #p.plot(band_chans, wnoise_rms[bl],':', label='$T_{wnoise}$')
    p.plot(band_chans, avg_trms, 'm--', label='$T_{avg,e-o}$')
    p.plot(band_chans, trms_the[bl],'--', label='$T_{RMS}$')
    p.plot(band_chans, trms_blavg, label='$T_{e-o,blavg}$')
    #p.plot(band_chans, noise_blavg,'-.', label='$T_{noise,blavg}$')
    #p.axhline(6.6408 *NOISE,color='k',linestyle=':')
    #p.plot(band_chans, wnoise_blavg,':', label='$T_{wnoise,blavg}$')
    p.plot(band_chans, trms_the_blavg,'--', label='$T_{RMS,blavg}$')
    p.plot(band_chans, theo_temp[bl],'k--')
    p.yscale('log')
    p.xlim([band_chans[0]-1,band_chans[-1]+1])
    p.legend(loc='upper right',ncol=3)
    p.title('%d_%d'%a.miriad.bl2ij(bl))
    #p.ylim([1e-2,1e4])
    if opts.plot:
        p.show()

    outfile = 'Trms_%d_%d.png'%a.miriad.bl2ij(bl)
    if not opts.output == '':
        outfile =opts.output+'/'+outfile
    p.savefig(outfile,format='png')
    p.clf()
