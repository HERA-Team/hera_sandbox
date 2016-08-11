#! /usr/bin/env python
import numpy as np, glob, sys
from matplotlib import pyplot as plt

JD = sys.argv[1] #this should be something like /data4/paper/2013EoR/2456678

#I wish you could dynamically name variables in Python
dxx1 = np.loadtxt(JD+'/fcal_xx_1/stats_1.txt')
dxx2 = np.loadtxt(JD+'/fcal_xx_2/stats_2.txt')

n1x,m1x,s1x = dxx1[:,0],dxx1[:,1],dxx1[:,2]
n2x,m2x,s2x = dxx2[:,0],dxx2[:,1],dxx2[:,2]

dyy1 = np.loadtxt(JD+'/fcal_yy_1/stats_1.txt')
dyy2 = np.loadtxt(JD+'/fcal_yy_2/stats_2.txt')

n1y,m1y,s1y = dyy1[:,0],dyy1[:,1],dyy1[:,2]
n2y,m2y,s2y = dyy2[:,0],dyy2[:,1],dyy2[:,2]

plt.errorbar(n1x,m1x,yerr=s1x,fmt='bo',ecolor='b',label='fc_xx_1')
plt.errorbar(n2x,m2x,yerr=s2x,fmt='co',ecolor='c',label='fc_xx_2')
plt.ylim(-100,100)
plt.xlim(-1,114)
plt.xlabel('Antenna #')
plt.ylabel('Delay')
plt.legend(loc='best')
plt.savefig(JD+'/stats_xx.png')
#plt.show()
plt.close()

plt.errorbar(n1y,m1y,yerr=s1y,fmt='ro',ecolor='r',label='fc_yy_1')
plt.errorbar(n2y,m2y,yerr=s2y,fmt='mo',ecolor='m',label='fc_yy_2')
plt.ylim(-100,100)
plt.xlim(-1,114)
plt.xlabel('Antenna #')
plt.ylabel('Delay')
plt.legend(loc='best')
plt.savefig(JD+'/stats_yy.png')
#plt.show()
plt.close()

npz_xx_1 = sorted(glob.glob(JD+'/fcal_xx_1/*npz'))
npz_xx_2 = sorted(glob.glob(JD+'/fcal_xx_2/*npz'))
npz_yy_1 = sorted(glob.glob(JD+'/fcal_yy_1/*npz'))
npz_yy_2 = sorted(glob.glob(JD+'/fcal_yy_2/*npz'))

Dx1,Dy1,Dx2,Dy2 = {},{},{},{}

Dlists = [Dx1,Dx2,Dy1,Dy2]
npzlists = [npz_xx_1,npz_xx_2,npz_yy_1,npz_yy_2]

for i in range(4):
    npz = np.load(npzlists[i][0])
    DD = Dlists[i]
    for k in npz.keys():
        if k.isdigit():
            DD[k] = []
    del(npz)
    npzs = npzlists[i]
    for npzfile in npzs:
        print '    Reading %s'%npzfile
        dn = np.load(npzfile)
        for k in dn.keys():
            if k.isdigit():
                try: DD[k].append(dn[k+'d'][0])
                except KeyError:
                    print 'KeyError for %s'%str(k)
                    continue

for k in Dx1.keys(): plt.plot(Dx1[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_xx_1',size=15)
plt.savefig(JD+'/fcal_xx_1/delay_per_file.png')
#plt.show()
plt.close()

for k in Dx2.keys(): plt.plot(Dx2[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_xx_2',size=15)
plt.savefig(JD+'/fcal_xx_2/delay_per_file.png')
#plt.show()
plt.close()

for k in Dy1.keys(): plt.plot(Dy1[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_yy_1',size=15)
plt.savefig(JD+'/fcal_yy_1/delay_per_file.png')
#plt.show()
plt.close()

for k in Dy2.keys(): plt.plot(Dy2[k])
plt.xlabel('File #')
plt.ylabel('Delay')
plt.suptitle('fc_yy_2',size=15)
plt.savefig(JD+'/fcal_yy_2/delay_per_file.png')
#plt.show()
plt.close()
