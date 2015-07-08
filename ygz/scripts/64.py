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

#c = 299792458.
#lamb = c/nu
#kB = 1.3806488E-23
#pref = (2*kB/lamb)*(2*kB/lamb)
pol,lst,ant,tauchan = opts.pol,opts.lst.split('_'),opts.ant.split('_'),int(opts.chan)
#dataDIR = '/Users/yunfanzhang/local/DATA64/'
dataDIR = ''
uv1 = a.miriad.UV(dataDIR+args[0])
uv2 = a.miriad.UV(dataDIR+args[1])

freqlist = n.array(uv1['sfreq'] + uv1['sdf']*n.arange(uv1['nchan']))
aa = a.cal.get_aa('psa6240_v003',freqlist)
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
pol = a.miriad.str2pol[pol]
t1,t2 = float(lst[0]),float(lst[1])
ant_dict, bl_dict = get_bls.get_bldict(aa)
bl_list1 = get_bls.get_equivalent(ant_dict,bl_dict,int(ant[0]),int(ant[1]))
bl_list2 = get_bls.get_equivalent(ant_dict,bl_dict,int(ant[2]),int(ant[3]))

#dt = uv1['inttime']
uv1.select('antennae',0,1,include=True)
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
    taulist = n.fft.fftfreq(int(uv1['nchan']),uv1['sdf'])
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
for item in sum:
    result[tau] = sum[tau]/cnt
    P.append(abs(sum[tau])/cnt)
print cnt,result[taulist[tauchan]]

plotp.P_v_Eta(taulist,P)
