__author__ = 'yunfanzhang'
import aipy as a, numpy as n, capo, pylab as p, capo
import  sys, export_beam, boot_simple
import delay_transform as dl_tr, plot_pspec as plotp
from capo import cosmo_units

c = a.const.c                 #cm/s
Mpc2m = 3.086E22 #meters
nu = n.arange(100, 200, 10)*1.E6
nu0 = 151.5
z, Omp, hub = 8.5, 0.63, 0.75
X2Y = capo.pspec.X2Y(z)      #[h^-3 Mpc^3] / [str * GHz]
bb = 0.1*(20./203)   # total window in GHz
df = 0.1/203

nchan = 20

DIR1 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/'
DIR2 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/'
Cname = 'P0.15.cue_mod'
uv1 = a.miriad.UV(DIR1+'pspec_2456249.20169.uv/')
uv2 = uv1
#uv2 = a.miriad.UV(DIR2+'pspec_0_38_2456249.20169.uv/')

freqflist = n.array((uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan'])))  #GHz
freqlist = n.array((uv1['sfreq'] + uv1['sdf']*92 + uv1['sdf']*n.arange(10)))  #GHz
aa = a.cal.get_aa('psa6240_v003',freqlist)
pol = a.miriad.str2pol['xx']
#taulist = n.fft.fftfreq(int(uv1['nchan']),uv1['sdf']*1000)
taulist = n.fft.fftfreq(nchan,uv1['sdf'])
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
with open(Cname, 'r') as f1:
    for line in f1:
        fnls = line.rstrip('\n').split(' ')
        if fnls[1]=='None' or fnls[2]=='None':
            print "No file for", fnls
            continue
        uv1 = a.miriad.UV(DIR1+str(fnls[1]))
        uv2 = a.miriad.UV(DIR2+str(fnls[2]))

        ant = line.split(' ')[0].split(',')[0].split('_')
        T1,T2 = line.split(' ')[0].split(',')[1].split('_')
        Ompp = float(line.split(' ')[3])                                          #Ompp computed for different separation in uv plane
        for i in range(4): ant[i] = int(ant[i])
        T1,T2 = float(T1),float(T2)
        #if ant[0] != 0 or ant[1] != 26 or ant[2] != 0 or ant[3] != 38: continue   #only simulated these baselines
        norm = X2Y*bb*Omp*Omp/Ompp    #1.E6 is K2 to mK2

        uv1.rewind()
        uv1.select('clear',0,0,include=True)
        uv1.select('antennae', ant[0], ant[1],include=True)
        uv1.select('polarization',pol,0,include=True)
        #uv1.select('chan',92,102,include=True)    #doesn't actually work
        uv1.select('time',T1-dt/1.9999,T1+dt/1.9999,include=True)
        #preamble, datnu = uv1.read()
        try: preamble, datnu = uv1.read()
        except(IOError):
            print "no data1 for", fnls, "interval", T1-dt/1.9999, T1+dt/1.9999
            continue
        #print 'bl1: ', preamble
        datnu = datnu[90:110]*norm           # only one of data1, data2 needs to be normalized such
        datatau = dl_tr.nu2tau(datnu)
        data1.append(n.array(datatau).transpose())
        #data1 = n.array(data1).transpose()


        uv2.rewind()
        uv2.select('clear',0,0,include=True)
        uv2.select('antennae',ant[2],ant[3],include=True)
        uv2.select('polarization',pol,0,include=True)
        #uv2.select('chan',92,102,include=True)
        uv2.select('time',T2-dt/1.9999,T2+dt/1.9999,include=True)
        try: preamble, datnu = uv2.read() # loop over time
        except(IOError):
            print "no data2 for", fnls, "interval", T2-dt/1.9999, T2+dt/1.9999
            continue
        datnu = datnu[90:110]
        #print "data lenth", len(datnu)
        datatau = dl_tr.nu2tau(datnu)
        data2.append(n.array(datatau).transpose())
        #print tauchan, datatau[tauchan]

data1, data2 = n.array(data1), n.array(data2)
print "shape of data1 (pt, tau):", data1.shape
print "shape of data2:", data2.shape
data = n.multiply(n.conjugate(data1), data2)
print "shape of data: ", data.shape
print data.all()

P=[]
for ind in n.arange(len(taulist)): P.append(n.mean(data.T[ind]))
#print len(taulist), len(P)
#kz = taulist*2*n.pi/Y
kz = cosmo_units.eta2kparr(taulist*1.E-9,z)     #This function needs tau in Hz^-1


#print "shapes of arrays:", data1.shape, data2.shape
#Bootstrap resampling
B = 100
bootmean, booterr = boot_simple.bootstrap(B, data)


#plotting
fig = p.figure()
ax = fig.add_subplot(311)
#plotp.P_v_Eta(ax,kz,P)
ax.set_xlabel('kz')
ax.set_ylabel(r'$P(k) K^{2} (h^{-1} Mpc)^{3}$')
p.plot(kz,P,'bo')
ax = fig.add_subplot(312)
ax.errorbar(kz, bootmean, yerr=booterr, fmt='ok', ecolor='gray', alpha=0.5)
#ax.set_ylim([0,0.5])
#ax.set_yscale('log')
ax.set_xlabel('kz')
ax.set_ylabel(r'$P(k) K^{2} (h^{-1} Mpc)^{3}$')
ax = fig.add_subplot(313)
plotp.Del_v_Eta(ax,kz,P)
p.show()

