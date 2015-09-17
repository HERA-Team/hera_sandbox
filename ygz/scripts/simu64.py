__author__ = 'yunfanzhang'
import aipy as a, numpy as n, capo, pylab as p, capo
import  sys, export_beam, boot_simple, os
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
schan = 90
polstr = 'xx'

DIR1 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/'
DIR2 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_38/'
F1, F2 = os.listdir(DIR1), os.listdir(DIR2)
for i in range(len(F1)): F1[i] = DIR1+F1[i]
for i in range(len(F2)): F2[i] = DIR2+F2[i]
Cname = 'P0.15.cue'
uv1 = a.miriad.UV(DIR1+'pspec_2456249.20169.uv/')
#uv2 = uv1
uv2 = a.miriad.UV(DIR2+'pspec_0_38_2456249.20169.uv/')

freqflist = n.array((uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan'])))  #GHz
freqlist = n.array((uv1['sfreq'] + uv1['sdf']*schan + uv1['sdf']*n.arange(nchan)))  #GHz
aa = a.cal.get_aa('psa6240_v003',freqflist)
aa.set_active_pol(polstr)
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15)
src.compute(aa)
pol = a.miriad.str2pol[polstr]
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

with open(Cname, 'r') as f1:
    #f1.read()    #skip header
    for line in f1:
        bl,DelT,Opp = line.rstrip('\n').split(',')
        if bl != '0_26_0_38' and bl != '0_38_0_26': continue
        DelT, Opp = float(DelT), float(Opp)
        norm = X2Y*bb*Omp*Omp/Opp    #1.E6 is K2 to mK2
        T1, dat1, flg1 = capo.arp.get_dict_of_uv_data(F1,antstr='0_26',polstr=polstr)
        T2, dat2, flg2 = capo.arp.get_dict_of_uv_data(F2,antstr='0_38',polstr=polstr)

        data1 = dat1[283][polstr]
        data2 = dat2[295][polstr]
        print T2[0], T1[0]

        while T1[0]+DelT>T2[0]:
            T2 = T2[1:]
            data2 = data2[1:]
        while len(T1)>len(T2):
            T1 = T1[:-1]
            data1 = data1[:-1]
        while T1[0]+DelT<T2[0]:
            T1 = T1[1:]
            data1 = data1[1:]
        while len(T2)>len(T1):
            T2 = T2[:-1]
            data2 = data2[:-1]

        print data1.shape, data2.shape  #(time,freq)

        phs1, phs2 = n.zeros(data1.shape,dtype='complex64'), n.zeros(data2.shape,dtype='complex64')
        t026 = 2456249.2666900107
        t038 = 2456249.3086900176

        ind = 0
        aa.set_jultime(t026)
        src.compute(aa)
        for t1 in T1:
            phs1[ind][:] = aa.gen_phs(src,0,26)
            ind = ind+1
        ind = 0
        aa.set_jultime(t038)
        src.compute(aa)
        for t2 in T2:
            phs2[ind][:] = aa.gen_phs(src,0,38)
            ind = ind+1
        data1, data2 = n.multiply(data1,phs1),n.multiply(data2,phs2)

        dau1, dau2 = [],[]
        for datnu in data1:
            datnu = datnu[schan:schan+nchan]*norm           # only one of data1, data2 needs to be normalized such
            datatau = dl_tr.nu2tau(datnu)
            dau1.append(datatau)
        for datnu in data2:
            datnu = datnu[schan:schan+nchan]*norm           # only one of data1, data2 needs to be normalized such
            datatau = dl_tr.nu2tau(datnu)
            dau2.append(datatau)


dau1, dau2 = n.array(dau1), n.array(dau2)
data = n.multiply(n.conjugate(dau1), dau2)
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

