#! /usr/bin/env python
import sys, numpy as np, aipy
from matplotlib import pyplot as plt
import os

uvfiles = sys.argv[1:]
#get initial info from first file
uv = aipy.miriad.UV(uvfiles[0])
nants,nchan = uv['nants'],uv['nchan']
uv.select('antennae',0,1,include=True)
_T=[]
for p,d in uv.all(): 
    _,_t,_ = p
    _T.append(_t)
ntimes = len(_T)*len(uvfiles)
del(uv)

if not os.path.exists('~/test.npz'):
    stor = np.zeros((nants,nchan,ntimes),dtype='complex128')

    for i in range(nants):
        print 'auto %i'%i
        counter=0
        for j,uvf in enumerate(uvfiles):
            uv = aipy.miriad.UV(uvf) #is there a way to undo a selection? would move this out of ant loop
            uv.select('antennae',i,i,include=True) #take auto
            for p,d in uv.all(): 
                stor[i,:,counter] = d
                counter+=1
            del(uv)
    np.savez('test.npz',stor=stor)    

print 'Plotting'
stor = np.load('test.npz')['stor']
#s = int(np.ceil(np.sqrt(nants)))
f,axarr = plt.subplots(12,11,sharex=True,sharey=True)
for i,ax in enumerate(axarr.ravel()):
    try:
        auto = np.abs(stor[i,:,:]) #np.log10(np.abs(stor[i,:,:]))
        line = np.nanmean(auto,axis=1)
        ax.plot(range(nchan),line)
        ax.text(170,25,str(i))
        #ax.set_ylim(0,35)
        #ax.set_xlim(0,203)
        #ax.legend()
    except IndexError:
        continue
        #ax.imshow(np.zeros_like(auto),cmap='binary')
    
plt.show()
