import numpy as np
import aipy as a
import matplotlib.pyplot as plt

uv = a.miriad.UV('/Users/JoshDillon/Desktop/PAPER/even/lst.2456242.20961.uvA')

generalInfo,data,flags = uv.read(raw=True)

print flags

# uv.rewind() goes back to the beginning

lst = []
flags = []
data = []
a.scripting.uv_selector(uv, '35_18', 'I') #baseline 1,4 and polarization I
for (crd,t,(i,j)),d,f in uv.all(raw=True):
    lst.append(t)
    flags.append(f)
    data.append(d)
    
plt.plot(lst,np.real(np.asarray(data)[:,30]))
plt.show()

a.plotuv.

#add to bashrc: add my scripts folder to my path
#airbears2-10-142-46-118:even JoshDillon$ export PATH=$PATH:/Users/JoshDillon/Desktop/capo/jsd/scripts
