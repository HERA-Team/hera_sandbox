__author__ = 'yunfanzhang'
import aipy as a, numpy as n, capo, pylab as p, capo
import optparse, sys, os, math
import get_bls, boot_simple
import delay_transform as dl_tr, plot_pspec as plotp
from capo import cosmo_units

#o = optparse.OptionParser()
#a.scripting.add_standard_options(o, chan=True, ant=True, pol=True)
#a.scripting.add_standard_options(o, ant=True, pol=True)
#o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time points to include')
#opts,args = o.parse_args(sys.argv[1:])
#print opts, args

c = a.const.c                 #cm/s
Mpc2m = 3.086E22 #meters
nu = n.arange(100, 200, 10)*1.E6
nu0 = 0.1515
z, Omp, hub, Ompp = 8.5, 0.63, 0.75, 0.31
X2Y = capo.pspec.X2Y(z)      #[h^-3 Mpc^3] / [str * GHz]
bb = 0.1*(20./203)   # total window in GHz
BB = bb*1.E9  #total window in Hz
df = 0.1/203

nchan = 20

#pol,lst,ant,tauchan = opts.pol,opts.lst.split('_'),opts.ant.split('_'),int(opts.chan)
pol = 'xx'
DIR = '/Users/yunfanzhang/local/DATA64/testDIR/'
uv1 = a.miriad.UV(DIR+'zen.2456249.24345.uvcRREcACOc')

freqflist = n.array((uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan'])))  #GHz
freqlist = n.array((uv1['sfreq'] + uv1['sdf']*92 + uv1['sdf']*n.arange(10)))  #GHz
aa = a.cal.get_aa('psa6240_v003',freqlist)
#src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
pol = a.miriad.str2pol[pol]

#taulist = n.fft.fftfreq(int(uv1['nchan']),uv1['sdf']*1000)
taulist = n.fft.fftfreq(nchan,uv1['sdf'])                                  #GHz
taulist = n.fft.ifftshift(taulist)
print uv1['sdf']

#dt = uv1['inttime']
uv1.select('antennae',0,26,include=True)
preamble, junk = uv1.read()
t01 = preamble[1]
uv1.read()
preamble, junk = uv1.read()
t02 = preamble[1]
dt = t02-t01


#bl_list1 = (0,26)]
#bl_list2 = [(0,26)]
bl1, bl2 = (0,26), (0,26)
ant_dict, bl_dict = get_bls.get_bldict(aa)
bl_list1 = get_bls.get_equivalent(ant_dict, bl_dict,0,26)

sum = []
data1,data2 = [],[]
cnt = 0
#norm = 1.E6*X2Y*bb*Omp*Omp/Ompp    #1.E6 is K to mK
norm = X2Y*bb*Omp*Omp/Ompp*capo.pspec.jy2T(nu0)**2   #jy2T returns mk/Jy for nu0 in GHz
for fn in os.listdir(DIR):
    print 'processing file', fn
    for bl1 in bl_list1:                                               #Each baseline is multiplied by itself
        bl2 = bl1
        uv1 = a.miriad.UV(DIR+fn)
        uv2 = uv1

        uv1.rewind()
        uv1.select('clear',0,0,include=True)
        uv1.select('antennae',bl1[0],bl1[1],include=True)
        uv1.select('polarization',pol,0,include=True)

        for preamble, datnu in uv1.all():
            print preamble
            datnu = datnu.filled(fill_value=0)[90:110]    #have to fill because fft ignores mask?
            datatau = dl_tr.nu2tau(datnu)
            data1.append(datatau.transpose())
        #data1 = n.array(data1).transpose()

        uv2.rewind()
        uv2.select('clear',0,0,include=True)
        uv2.select('antennae',bl2[0],bl2[1],include=True)
        uv2.select('polarization',pol,0,include=True)

        for preamble, datnu in uv2.all():  # loop over time
            datnu = datnu.filled(fill_value=0)[90:110]
            datatau = dl_tr.nu2tau(datnu)
            data2.append(datatau.transpose())
            #print tauchan, datatau[tauchan]


#print "data shapes", data1.shape, data2.shape
print "Average over %d time points" % len(data1)
data1, data2 = n.ma.array(data1), n.ma.array(data2)
data = n.multiply(n.conjugate(data1), data2)*norm
#data is shaped [timesample, channel]
P = n.ma.array([])
for ind in n.arange(len(taulist)): P = n.append(P,n.mean(data.T[ind]))
print P
#print len(taulist), len(P)
#kz = taulist*2*n.pi/Y
kz = cosmo_units.eta2kparr(taulist*1.E-9,z)     #This function needs tau in Hz^-1


#print "shapes of arrays:", data1.shape, data2.shape
#Bootstrap resampling
B = 10
bootmean, booterr = boot_simple.bootstrap(B, data)


#plotting
fig = p.figure()
ax = fig.add_subplot(311)
#plotp.P_v_Eta(ax,kz,P)
ax.set_xlabel('kz')
ax.set_ylabel(r'$P(k) mK^{2} (h^{-1} Mpc)^{3}$')
p.plot(kz,P,'bo')
ax.set_yscale('log')
ax = fig.add_subplot(312)
ax.errorbar(kz, bootmean, yerr=booterr, fmt='ok', ecolor='gray', alpha=0.5)
#ax.set_ylim([0,0.5])
ax.set_yscale('log')
ax.set_xlabel('kz')
ax.set_ylabel(r'$P(k) mK^{2} (h^{-1} Mpc)^{3}$')
ax = fig.add_subplot(313)
plotp.Del_v_Eta(ax,kz,P)
p.show()