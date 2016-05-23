import numpy as n,sys,os
from pylab import *

maxsubplots = 6*5
for maskfile in sys.argv[1:]:
    mask = n.load(maskfile)
    figure()
    maxsubplots = min((len(mask.files),maxsubplots))
    rows = n.ceil(n.sqrt(maxsubplots))
    cols = n.ceil(maxsubplots/rows)
    plotnum = 1
    for ant in mask:
        if plotnum>maxsubplots: 
            plotnum=1
            show()
        subplot(rows,cols,plotnum)
        imshow(mask[ant],interpolation='nearest',aspect='auto')
        title(ant)
        plotnum += 1
        print n.sum(mask[ant]-mask['0x'])
    suptitle(maskfile)
    show()
