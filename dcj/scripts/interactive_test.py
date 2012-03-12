#!/usr/bin/env python
#
#  interactive_test.py
#  
#
#  Created by Danny Jacobs on 9/3/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#
from pylab import *
import time

ion()

tstart = time.time()               # for profiling
x = arange(0,2*pi,0.01)            # x-array
line, = plot(x,sin(x))
for i in arange(1,200):
    line.set_ydata(sin(x+i/10.0))  # update the data
    draw()                         # redraw the canvas

print 'FPS:' , 200/(time.time()-tstart)