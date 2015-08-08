__author__ = 'yunfanzhang'
import aipy as a, numpy as n, capo, pylab as p, capo
import optparse, sys, random, os
import get_bls
import delay_transform as dl_tr, plot_pspec as plotp

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True, ant=True, pol=True)
#a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time points to include')
opts,args = o.parse_args(sys.argv[1:])
print opts, args

c = 299792458.
Mpc2m = 3.086E22 #meters
nu = n.arange(100, 200, 10)*1.E6
nu0 = 150
z, Omp, hub, Ompp = 8.5, 0.63, 0.75, 0.31
X2Y = capo.pspec.X2Y(z)      #[h^-3 Mpc^3] / [str * GHz]
bb = 0.1*(20./203)   # total window in GHz
df = 0.1/203

nchan = 20

pol,lst,ant,tauchan = opts.pol,opts.lst.split('_'),opts.ant.split('_'),int(opts.chan)
DIR = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/'
dataDIR = ''
uv1 = a.miriad.UV(dataDIR+args[0])
uv2 = a.miriad.UV(dataDIR+args[1])

freqflist = n.array((uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan']))*1000)  #MHz
freqlist = n.array((uv1['sfreq'] + uv1['sdf']*92 + uv1['sdf']*n.arange(10))*1000)  #MHz
aa = a.cal.get_aa('psa6240_v003',freqlist)
#src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
pol = a.miriad.str2pol[pol]
t1,t2 = float(lst[0]),float(lst[1])
#bl_list1 = (0,26)]
#bl_list2 = [(0,26)]
bl1, bl2 = (0,26), (0,26)
#taulist = n.fft.fftfreq(int(uv1['nchan']),uv1['sdf']*1000)
taulist = n.fft.fftfreq(nchan,uv1['sdf']*1000)
taulist = n.fft.ifftshift(taulist)

#dt = uv1['inttime']
uv1.select('antennae',0,26,include=True)
preamble, junk = uv1.read()
t01 = preamble[1]
uv1.read()
preamble, junk = uv1.read()
t02 = preamble[1]
dt = t02-t01

sum = []
data1,data2 = [],[]
cnt = 0
norm = 1.E6*X2Y*bb*Omp*Omp/Ompp    #1.E6 is K to mK
for fn in os.listdir(DIR):
    uv1 = a.miriad.UV(DIR+fn)
    uv2 = a.miriad.UV(DIR+fn)

    uv1.rewind()
    uv1.select('clear',0,0,include=True)
    uv1.select('antennae',bl1[0],bl1[1],include=True)
    uv1.select('polarization',pol,0,include=True)
    #uv1.select('chan',92,102,include=True)    #doesn't actually work
    #uv1.select('time',t1-dt/2.001,t1+dt/2.001,include=True)
    #preamble, datnu = uv1.read()
    for preamble, datnu in uv1.all():
        print 'bl1: ', preamble
        datnu = datnu[90:110]*norm
        datatau = dl_tr.nu2tau(datnu)
        data1.append(n.array(datatau).transpose())
    #data1 = n.array(data1).transpose()

    uv2.rewind()
    uv2.select('clear',0,0,include=True)
    uv2.select('antennae',bl2[0],bl2[1],include=True)
    uv2.select('polarization',pol,0,include=True)
    #uv2.select('chan',92,102,include=True)
    #uv2.select('time',t2-dt/2.001,t2+dt/2.001,include=True)
    #preamble, datnu = uv2.read() # loop over time
    for preamble, datnu in uv2.all():  # loop over time
        datnu = datnu[90:110]
        print "data lenth", len(datnu)
        datatau = dl_tr.nu2tau(datnu)
        data2.append(n.array(datatau).transpose())
        #print tauchan, datatau[tauchan]


#print "data shapes", data1.shape, data2.shape
print "Average over %d time points" % len(data1)
for ind in range(len(taulist)): sum.append(0.)
for ine in range(len(data1)):
    for ind in range(len(taulist)):
        tau = taulist[ind]
        sum[ind] = sum[ind] + n.sum(n.conjugate(data1[ine][ind])*data2[ine][ind])
    cnt = cnt + len(data1[ine])

result = {}
P=[]
for ind in n.arange(len(sum)):
    #result[taulist[ind]] = sum[ind]/cnt
    #P.append(abs(sum[ind])/cnt/pref*XSY/B/W*1.E-52*1.E12)   #1Jy=E-26W/m2/Hz
    P.append(sum[ind]/cnt)
print len(taulist), len(P)
#kz = taulist*2*n.pi/Y
kz = capo.cosmo_units.eta2kparr(taulist,z)
plotp.P_v_Eta(kz,P)


#Bootstrap resampling

B = 1000
data1 = n.array(data1)
print data1.shape
boot = []
for b in range(B):
    temps = []
    for i in range(len(data1[0])):    #for each channel
        temps.append(0.)
        for j in range(len(data1)):      #for each sample
            poin = random.choice(data1)[i]    #for each  b, channel i, resample from data1
            temps[i] = temps[i] + n.conjugate(poin)*poin/len(data1)    #temp is the b^th PS estimate
    boot.append(temps)
boot = n.array(boot).transpose()
print boot.shape
bootmean,booterr = [],[]
for ch in range(len(boot)):
    mean = n.sum(boot[ch])/B
    sig = n.sqrt(n.sum((boot[ch]-mean)**2)/(B-1))
    bootmean.append(mean)
    booterr.append(2*sig)
#p.hist(boot[5])
#p.show()

fig, ax = p.subplots()
ax.errorbar(kz, bootmean, yerr=booterr, fmt='ok', ecolor='gray', alpha=0.5)
ax.set_ylim([0,0.5])
#ax.set_yscale('log')
ax.set_xlabel('kz')
ax.set_ylabel(r'$P(k) K^{2} (h^{-1} Mpc)^{3}$')
p.show()