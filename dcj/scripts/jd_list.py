#!/bin/env python
import string,sys
filestring=sys.argv[1]
A=[string.join(s.split('.')[1:3],'.') for s in sys.argv[1:]]
for a in A:print a,' ',

