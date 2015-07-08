#! /usr/bin/env python
import aipy as a
import sys
VERB=False
for f in sys.argv[1:]:
    print f,
    sys.stdout.flush()
    try:
        UV = a.miriad.UV(f)
        for all in UV.all():
            continue
        print "[OK]"
        sys.stdout.flush()
    except:
        print "[ERROR]"
        if VERB: 
            print sys.exc_info()[0]
            sys.stdout.flush()

