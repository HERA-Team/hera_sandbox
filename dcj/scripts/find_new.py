#!/usr/bin/env python
#
#  find_new.py
#  
#
#  Created by Danny Jacobs on 6/9/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import glob

uvcbs = set(['.'.join(s.split('.')[1:3]) for s in glob.glob('*.uvcb')])
zips = set(['.'.join(s.split('.')[1:3]) for s in glob.glob('*.bz2')])
new = zips.difference(uvcbs)
for i in range(len(new)): print('zen.%s.uv.tar.bz2'%(new.pop(),))
