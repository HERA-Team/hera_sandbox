#!/usr/bin/env python
import sys

files = sys.argv[1:]
files1 = files[1::3]
files2 = files[2::3]
files3 = files[0::3]

print files1,files2,files3
