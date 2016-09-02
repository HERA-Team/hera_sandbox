import numpy as np, glob, sys, os
from matplotlib import pyplot as plt

def get_results(JDs,rnd=1,epoch_bounds=[2456674]):
    results = np.zeros((112,len(JDs),2,2)) #antnum, day, mean/sd, epoch 1/2
    ep = 0
    for i,dir in enumerate(JDs):
        print '%s'%dir,i
        statfile = dir+'/fcal_xx_%i/stats_%i.txt'%(rnd,rnd)
        if not os.path.exists(statfile):
            print '%s does not exist'%statfile
            continue
        #if int(dir.split('/')[-1]) in epoch_bounds: ep+=1
        N,m,s = np.loadtxt(statfile,unpack=True)
        N = map(int,N)
        for j,antnum in enumerate(N):
            results[antnum,i,0,ep] = m[j]
            results[antnum,i,1,ep] = s[j]
    return results

JDs = sys.argv[1:]
r1,r2 = get_results(JDs,rnd=1),get_results(JDs,rnd=2)

plt.imshow(np.log10(r1[:,:,1,0]),aspect='auto',interpolation='None')
plt.xlabel('Day #')
plt.ylabel('Antenna #')
plt.suptitle('fcal_xx_1')
plt.colorbar()
plt.show()

plt.imshow(np.log10(r2[:,:,1,0]),aspect='auto',interpolation='None')
plt.xlabel('Day #')
plt.ylabel('Antenna #')
plt.suptitle('fcal_xx_2')
plt.colorbar()
plt.show()


