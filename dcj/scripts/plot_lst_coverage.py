#! /usr/bin/env python
from pylab import *
import numpy as n,sys
m=1
for file in sys.argv[1:]:
    D = n.loadtxt(file)
    name=file.split('_')[0]
    if name.startswith('total'):plot(D[:,0],D[:,1]*10,'k',lw=3,label=name)
    else:plot(D[:,0],D[:,1]*10,label=name)
    if m<D[:,1].max()*10: m=D[:,1].max()*10
legend(loc='best')
xlabel('local sidereal time [h]')
ylabel('time [min]')
ylim([0,m*1.2])
show()

