#!/usr/bin/env python
import numpy as np, matplotlib.pyplot as plt, sys

npzs = sys.argv[1:]
print npzs
storage = np.zeros((112,len(npzs)*19,203),dtype='complex128')
f,axarr = plt.subplots(8,14)
axs = axarr.ravel()

for i,npz in enumerate(npzs):
    data = np.load(npz)   
    print 'Working on %s'%npz     
    for ant in range(0,112):
        try: storage[ant,i*19:(i+1)*19,:] = data[str(ant)+'x'] #XXX limited to xx pol fc run
        except KeyError: 
            print 'Key Error for antenna %i in %s'%(ant,npz)
            storage[ant,i*19:(i+1)*19,:] = np.zeros((19,203),dtype='complex128')

for i,ax in enumerate(axs):
    ax.imshow(np.angle(storage[i,:,:]))
    ax.set_title(str(i)+'x')
plt.show()
    
