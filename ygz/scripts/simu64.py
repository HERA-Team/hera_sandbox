__author__ = 'yunfanzhang'
import aipy as a, numpy as n, capo
import optparse, sys
import get_bls
import delay_transform as dl_tr, plot_pspec as plotp

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True, ant=True, pol=True)
#a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time points to include')
opts,args = o.parse_args(sys.argv[1:])
print opts, args

c = 299792458.
nu0 = 150.E6  #the mean frequency 150 MHz
lamb = c/nu0
kB = 1.3806488E-23
pref = (2*kB/lamb*lamb)*(2*kB/lamb*lamb)
z, Omm, hub = 8.5, 0.27, 0.75
#Y = 17 (((1+z)/10)/(Omm*hub*hub/0.15))^0.5 #Mpc/MHz
#X = 1.9 ((1+z)/10)^0.2/hub                 # Mpc/arcmin
Y = 17*(((1+z)/10)/(Omm/0.15))**0.5 #h Mpc/MHz
X = 1.9*((1+z)/10)**0.2                # h Mpc/arcmin
B = 100   #MHz
W = 0.31


pol,lst,ant,tauchan = opts.pol,opts.lst.split('_'),opts.ant.split('_'),int(opts.chan)
#dataDIR = '/Users/yunfanzhang/local/DATA64/'
dataDIR = ''
uv1 = a.miriad.UV(dataDIR+args[0])
uv2 = a.miriad.UV(dataDIR+args[1])

freqlist = n.array((uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan']))*1000)  #MHz
aa = a.cal.get_aa('psa6240_v003',freqlist)
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
pol = a.miriad.str2pol[pol]
t1,t2 = float(lst[0]),float(lst[1])
bl_list1 = [(0,26)]
bl_list2 = [(0,26)]
taulist = n.fft.fftfreq(int(uv1['nchan']),uv1['sdf']*1000)
taulist = n.fft.fftshift(taulist)

#dt = uv1['inttime']
uv1.select('antennae',0,26,include=True)
preamble, junk = uv1.read()
t01 = preamble[1]
uv1.read()
preamble, junk = uv1.read()
t02 = preamble[1]
dt = t02-t01

cnt = 0
sum = []
data1,data2 = [],[]

for bl1 in bl_list1:
    data1 = []
    uv1.rewind()
    uv1.select('clear',0,0,include=True)
    uv1.select('antennae',bl1[0],bl1[1],include=True)
    uv1.select('polarization',pol,0,include=True)
    uv1.select('time',t1-dt/2.001,t1+dt/2.001,include=True)
    try: preamble, datnu = uv1.read()
    except(IOError): continue
    print 'bl1: ', preamble
    datatau = dl_tr.nu2tau(datnu)
    data1.append(datatau)
    print tauchan, datatau[tauchan]
    data1 = n.array(data1).transpose()

    for bl2 in bl_list2:
        data2 = []
        uv2.rewind()
        uv2.select('clear',0,0,include=True)
        uv2.select('antennae',bl2[0],bl2[1],include=True)
        uv2.select('polarization',pol,0,include=True)
        uv2.select('time',t2-dt/2.001,t2+dt/2.001,include=True)
        try: preamble, datnu = uv2.read() # should be a trivial loop, onlye has one entry
        except(IOError):
            print "Error reading ant (%d, %d), time (%f,%f)" % (bl2[0],bl2[1],t2-dt/2.001,t2+dt/2.001)
            continue
        #print 'bl2:', preamble
        print "data lenth", len(datnu)
        datatau = dl_tr.nu2tau(datnu)
        data2.append(datatau)
        #print tauchan, datatau[tauchan]

        data2 = n.array(data2).transpose()
        #print "data shapes", data1.shape, data2.shape
        for ind in range(len(taulist)):
            tau = taulist[ind]
            sum.append(n.sum(n.conjugate(data1[ind])*data2[ind]))
        cnt = cnt + len(data1)

result = {}
P=[]
for ind in n.arange(len(sum)):
    result[taulist[ind]] = sum[ind]/cnt
    P.append(abs(sum[ind])/cnt/pref*X*X*Y/B/W)
print cnt,result[taulist[tauchan]]
print len(taulist), len(P)
kz = taulist*2*n.pi/Y
plotp.P_v_Eta(kz,P)
