#!/usr/bin/env python
#
#  find_difference.py
#  
#
#  Created by Danny Jacobs on 6/9/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import glob,sys
"""
find_difference.py <type_A> <type_B>
Compare JD on the filenames in current directory for two kinds of files. Outputs 
files of type B for where there are no files of type A. Types can be any regular
expression as valid input to the python module glob.
ex:
find_difference.py '*.uvcb' '*.bz2'  #outputs the bz2 files for which there are no 
                                 #uvcb's 
"""

if len(sys.argv[1:])!=2: raise SyntaxError('Please input exactly two file types')
typeA = sys.argv[1]
typeB = sys.argv[2]
#print typeA,typeB
pref = glob.glob(typeA)[0].split('.')[0]
uvcbs = set(['.'.join(s.split('.')[1:3]) for s in glob.glob(typeA)])
zips = set(['.'.join(s.split('.')[1:3]) for s in glob.glob(typeB)])
new = zips.difference(uvcbs)
for i in range(len(new)): print('%s.%s*%s'%(pref,new.pop(),typeB))
