import numpy as n
from pylab import *
import sys,os,numpy as n
from capo.dcj import file2jd

chisq=[]
nprefix=3
media = []
means = []
jds = [file2jd(filename) for filename in sys.argv[1:]]
I = n.argsort(jds)
for datafile in [sys.argv[1:][i] for i in I]:
    data = np.fromfile(datafile, dtype='float32')
    nf = int(data[nprefix - 1])
    nt = len(data) / (nf + nprefix)
    if len(data) != nt * (nf + nprefix):
        raise IOError("File %s is not a valid omnichisq file."%datafile)
    
    data.shape = (nt, (nf + nprefix))
    data = data[:, nprefix:]
    chisq.append(data)
    media.append(n.median(data,axis=0))
    means.append(n.mean(data,axis=0))
#    print '.',
    print datafile
chisq = n.concatenate(chisq)
media = n.vstack(media)
means = n.vstack(means)
print chisq.shape,media.shape
imshow(n.log(chisq),aspect='auto',vmin=-6,vmax=6)
colorbar()
#figure()
#imshow(n.log(media),aspect='auto',vmin=-6,interpolation='nearest',vmax=-3)
#title('chisq median per file')
#colorbar()
#figure()
#imshow(n.log(means),aspect='auto',vmin=-6,interpolation='nearest',vmax=-3)
#colorbar()
#figure()
#for i in range(100):
#    plot(n.log(means[i,:])+n.std(means)*i)
show()
