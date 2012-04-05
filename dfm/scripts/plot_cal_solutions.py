import numpy as np
from pylab import *

calfile = 'calSTDs.npz'
nocalfile = 'NOcalSTDs.npz'

cal = np.load(calfile)
nocal = np.load(nocalfile)

figure(0)

subplot(231)
imshow(nocal['amp'],aspect='auto',interpolation='nearest',vmin=0,vmax=0.5)
colorbar()
title('No Calibration')
ylabel('AMP')

subplot(232)
imshow(cal['amp'],aspect='auto',interpolation='nearest',vmin=0,vmax=0.5)
colorbar()
title('Calibration')

subplot(233)
imshow(nocal['amp']-cal['amp'],aspect='auto',interpolation='nearest',vmin=-0.2,vmax=0.2)
colorbar()
title('NoCal - Cal')

subplot(234)
imshow(nocal['phs'],aspect='auto',interpolation='nearest',vmin=0,vmax=2)
colorbar()
ylabel('PHS')

subplot(235)
imshow(cal['phs'],aspect='auto',interpolation='nearest',vmin=0,vmax=2)
colorbar()

subplot(236)
imshow(nocal['phs']-cal['phs'],aspect='auto',interpolation='nearest',vmin=-0.2,vmax=2)
colorbar()

draw()

show()
