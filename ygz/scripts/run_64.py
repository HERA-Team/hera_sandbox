__author__ = 'yunfanzhang'
import aipy as a, numpy as n
import optparse, sys
import get_bls

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True, ant=True, pol=True)
o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time points to include')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa('psa6240_v003',n.array([.150739]))
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
pol = a.miriad.str2pol(opts.pol)
t1,t2 = opts.lst[0],opts.lst[1]
S = opts.chan
ant_dict, bl_dict = get_bls.get_bldict(aa)
bl_list1 = get_bls.get_equivalent(ant_dict,bl_dict,opts.ant[0][0],opts.ant[0][1])
bl_list2 = get_bls.get_equivalent(ant_dict,bl_dict,opts.ant[1][0],opts.ant[1][1])

print opts.ant, opts.lst
uv1 = a.miriad.UV(args[0])
uv2 = a.miriad.UV(args[1])

uv1.select('antennae',0,1,include=True)
preamble, junk = uv1.read()
t01 = preamble[1]
preamble, junk = uv1.read()
t02 = preamble[1]
dt = t02-t01
uv1.rewind()

cnt = 0
sum = 0
for bl1 in bl_list1:
    for bl2 in bl_list2:
        uv1.select('clear',0,0,include=True)
        uv1.select('antennae',bl1[0],bl1[1],include=True)
        uv1.select('polarization',pol,0,include=True)
        uv1.select('time',t1-dt/4.,t1+dt/4.,include=True)
        preamble, data1 = uv1.all()

        uv2.select('clear',0,0,include=True)
        uv2.select('antennae',bl2[0],bl2[1],include=True)
        uv2.select('polarization',pol,0,include=True)
        uv2.select('time',t2-dt/4.,t2+dt/4.,include=True)
        preamble, data2 = uv2.all()

        sum = sum + n.conjugate(data1[S])*data2[S]
        cnt = cnt + 1

result = sum/cnt
print cnt,result

