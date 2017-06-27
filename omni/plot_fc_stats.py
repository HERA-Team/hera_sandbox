#! /usr/bin/env python

import numpy as np, glob, sys
from matplotlib import pyplot as plt
import optparse

o = optparse.OptionParser()
o.set_usage('plot_fc_stats.py [options] /path/to/JD')
o.add_option('--suffix',dest='suffix',type='string',default='',help='This script will assume a psa128-type directory structure, with /path/to/JD/*uvcRRE files and /path/to/JD/fcal_xx_N/stats_N.txt files. A suffix will land after the fcal_xx_N, e.g. suffix="_allEW".')
o.add_option('--round3',dest='round3',action='store_true',help='If there was a third round of firstcal-ing, choose this. Otherwise it assumes only 2 rounds.')

opts,args = o.parse_args(sys.argv[1:])
assert(len(args)==1)
JD = args[0] #this should be something like /data4/paper/2013EoR/2456678

#stoopid python
fcaldirs=['fcal_xx_1','fcal_xx_2','fcal_yy_1','fcal_yy_2']
if opts.round3: fcaldirs.extend(['fcal_xx_3','fcal_yy_3'])
if len(opts.suffix)>0:
    _fcaldirs = []
    for x in fcaldirs:
        _fcaldirs.append(x+opts.suffix)
    fcaldirs = _fcaldirs

#I wish you could dynamically name variables
# ARP: you can: eval

dxx1 = np.loadtxt(JD+'/%s/stats_1.txt'%fcaldirs[0])
dxx2 = np.loadtxt(JD+'/%s/stats_2.txt'%fcaldirs[1])
dyy1 = np.loadtxt(JD+'/%s/stats_1.txt'%fcaldirs[2])
dyy2 = np.loadtxt(JD+'/%s/stats_2.txt'%fcaldirs[3])
if opts.round3:
    dxx3 = np.loadtxt(JD+'/%s/stats_3.txt'%fcaldirs[4])
    dyy3 = np.loadtxt(JD+'/%s/stats_3.txt'%fcaldirs[5])
    

n1x,m1x,s1x = dxx1[:,0],dxx1[:,1],dxx1[:,2]
n2x,m2x,s2x = dxx2[:,0],dxx2[:,1],dxx2[:,2]
n1y,m1y,s1y = dyy1[:,0],dyy1[:,1],dyy1[:,2]
n2y,m2y,s2y = dyy2[:,0],dyy2[:,1],dyy2[:,2]
if opts.round3:
    n3x,m3x,s3x = dxx3[:,0],dxx3[:,1],dxx3[:,2]
    n3y,m3y,s3y = dyy3[:,0],dyy3[:,1],dyy3[:,2]

plt.errorbar(n1x,m1x,yerr=s1x,fmt='bo',ecolor='b',label='fc_xx_1')
plt.errorbar(n2x,m2x,yerr=s2x,fmt='co',ecolor='c',label='fc_xx_2')
plt.errorbar(n3x,m3x,yerr=s3x,fmt='go',ecolor='g',label='fc_xx_3')

plt.ylim(-100,100)
plt.xlim(-1,114)
plt.xlabel('Antenna #')
plt.ylabel('Delay')
plt.legend(loc='best')
#plt.savefig(JD+'/stats_xx.png')
plt.show()
#plt.close()

plt.errorbar(n1y,m1y,yerr=s1y,fmt='ro',ecolor='r',label='fc_yy_1')
plt.errorbar(n2y,m2y,yerr=s2y,fmt='mo',ecolor='m',label='fc_yy_2')
plt.errorbar(n3y,m3y,yerr=s3y,fmt='yo',ecolor='y',label='fc_yy_3')

plt.ylim(-100,100)
plt.xlim(-1,114)
plt.xlabel('Antenna #')
plt.ylabel('Delay')
plt.legend(loc='best')
#plt.savefig(JD+'/stats_yy.png')
plt.show()
#plt.close()

#npz_xx_1 = sorted(glob.glob(JD+'/fcal_xx_1/*npz'))
#npz_xx_2 = sorted(glob.glob(JD+'/fcal_xx_2/*npz'))
#npz_yy_1 = sorted(glob.glob(JD+'/fcal_yy_1/*npz'))
#npz_yy_2 = sorted(glob.glob(JD+'/fcal_yy_2/*npz'))

npz_xx_1 = sorted(glob.glob(JD+'/*npz'))
npz_xx_2 = sorted(glob.glob(JD+'/*npz'))
npz_yy_1 = sorted(glob.glob(JD+'/*npz'))
npz_yy_2 = sorted(glob.glob(JD+'/*npz'))


npz_xx_1 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[0]))
npz_xx_2 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[1]))
npz_yy_1 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[2]))
npz_yy_2 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[3]))
if opts.round3:
    npz_xx_3 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[4]))
    npz_yy_3 = sorted(glob.glob(JD+'/%s/*npz'%fcaldirs[5]))

Dx1,Dy1,Dx2,Dy2,Dx3,Dy3 = {},{},{},{},{},{}

Dlists = [Dx1,Dx2,Dy1,Dy2]
if opts.round3: Dlists.extend([Dx3,Dy3])
npzlists = [npz_xx_1,npz_xx_2,npz_yy_1,npz_yy_2]
if opts.round3: npzlists.extend([npz_xx_3,npz_yy_3])

if not opts.round3: NN=4
else: NN=6

for i in range(NN):
    npz = np.load(npzlists[i][0])
    DD = Dlists[i]
    for k in npz.keys():
#        if k.isdigit():
        if k.startswith('d'):
            k=k[1:]
#            print 'loaded ',k
            DD[k] = []
    del(npz)
    npzs = npzlists[i]
    for npzfile in npzs:
        print '    Reading %s'%npzfile
        dn = np.load(npzfile)
        for k in dn.keys():
#            if k.isdigit():
            if k.startswith('d'):
                k=k[1:]
#                print 'bfnadsfknadslkfn ',k
                try: DD[k].append(dn['d'+k])
                except KeyError:
                    print 'KeyError for %s'%str(k)
                    continue

for k in Dx1.keys(): plt.plot(Dx1[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_xx_1',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[0])
#plt.show()
plt.close()

for k in Dx2.keys(): plt.plot(Dx2[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_xx_2',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[1])
#plt.show()
plt.close()

for k in Dx3.keys(): plt.plot(Dx3[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_xx_3',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[4])
#plt.show()
plt.close()

for k in Dy1.keys(): plt.plot(Dy1[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_yy_1',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[2])
#plt.show()
plt.close()

for k in Dy2.keys(): plt.plot(Dy2[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_yy_2',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[3])
#plt.show()
plt.close()

for k in Dy3.keys(): plt.plot(Dy3[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_yy_3',size=15)
plt.savefig(JD+'/%s/delay_per_file.png'%fcaldirs[5])
#plt.show()
plt.close()
