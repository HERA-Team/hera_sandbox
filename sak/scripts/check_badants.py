#! /usr/bin/env python
import numpy as np, glob, sys
from matplotlib import pyplot as plt

files = sys.argv[1:]
spl = files[0].split('_')
p = spl[1]
N = spl[2][0]
master = []

checkerboard = np.zeros((len(files),111))

for day,fname in enumerate(files):
    with open(fname) as f:
        content = f.readlines()
        l=content[0].split(' ')
        JD,b = int(l[0]),l[1].split(',')
        if day<43 and day > 39: print JD
        badants = [x.strip('\n') for x in b]
        badants = map(int, badants)
        for ba in badants: checkerboard[day,ba]=1
        master+=badants
m = np.array(master)
count = np.bincount(m)
x = range(len(count))
print count
fig,ax = plt.subplots()

ax.step(x,count,where='mid')
#ax.fill_between(x,0,count,step_where='mid',facecolor='blue',alpha=0.5)
ax.set_xlim(0,111)
ax.set_ylim(0,len(files))
ax.set_xlabel('Antenna Number')
ax.set_ylabel('Badness count over %i days'%len(files))
ax.set_title('%s round %s'%(p,N))
#plt.show()
plt.close()

plt.imshow(checkerboard,aspect='auto',interpolation='None',cmap='binary')
plt.xlabel('Antenna number')
plt.ylabel('Night number')
plt.suptitle('%s round %s'%(p,N))
plt.show()


np.savez('badants_occ_%s_%s.npz'%(p,N),occ=count)



        
