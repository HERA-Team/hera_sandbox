#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import atpy

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

f0 = n.load(args[0])
inpnums = f0['inpnums'].flatten()
del(f0)

print args

C,C_gain,times = {},{},[]
for fname in args:
    f = n.load(fname)
    times.append(f['time'])
    for ind, inp in enumerate(inpnums):
        if not C.has_key(inp):
            C[inp] = []
            C_gain[inp] = []
        C[inp].append(f['C'].flatten()[ind])
        C_gain[inp].append(f['C_gain'].flatten()[ind])

times = n.array(times)
aa = a.cal.get_aa(opts.cal,n.array([.150]))
lst = []
cat = a.src.get_catalog(srcs=['Sun','sgr'])
alt_sun,az_sun = [],[]
alt_sgr,az_sgr = [],[]
for t in times:
    aa.set_jultime(t)
    cat.compute(aa)
    lst.append(aa.sidereal_time() * (12/n.pi))
    alt_sun.append(cat['Sun'].alt * (180/n.pi))
    az_sun.append(cat['Sun'].az * (180/n.pi))
    alt_sgr.append(cat['sgr'].alt * (180/n.pi))
    az_sgr.append(cat['sgr'].az * (180/n.pi))

if False:
    for inp in inpnums:
        p.subplot(121)
        p.title('Delay (ns)')
        p.plot(times,C[inp],'.',label=str(inp))
        p.subplot(122)
        p.title('Gain (dB)')
        p.plot(times,10*n.array(C_gain[inp]),'.',label=str(inp))
    #p.legend()
    p.show()

filename = 'cal3_prms.vot'
if True:
    t = atpy.Table()
    ants = C.keys()
    ants.sort()
    for ant in ants:
        t.add_column('A%d_Dly' % ant, C[ant], dtype='<f8')
        t.add_column('A%d_Gain' % ant, C_gain[ant], dtype='<f8')
    t.add_column('Time', times, dtype='<f8', unit='JD')
    t.add_column('LST', lst, dtype='<f8', unit='hours')
    t.add_column('Alt_Sun', alt_sun, dtype='<f8',unit='deg')
    t.add_column('Az_Sun', az_sun, dtype='<f8',unit='deg')
    t.add_column('Alt_sgr', alt_sgr, dtype='<f8',unit='deg')
    t.add_column('Az_sgr', az_sgr, dtype='<f8',unit='deg')
    t.write(filename,overwrite=True)

if False:
    Cs = n.zeros((len(C.keys())+1,len(C.values()[0])+1))
    Cgs = n.zeros((len(C_gain.keys())+1,len(C_gain.values()[0])+1))
    keys = C.keys()
    keys.sort()
    Cs[:,0] = n.append(n.array([0]),keys)
    Cgs[:,0] = n.append(n.array([0]),keys)
    Cs[0,:] = n.append(n.array([0]),times)
    Cgs[0,:] = n.append(n.array([0]),times)
    for i,j in enumerate(keys):
        Cs[i+1,1:] = n.around(C[j],2)
        Cgs[i+1,1:] = n.around(10*n.array(C_gain[j]),2)

    n.savetxt('C3_dly_dat.txt',Cs.T,fmt='%5.4f')
    n.savetxt('C3_gain_dat.txt',Cgs.T,fmt='%5.4f')
