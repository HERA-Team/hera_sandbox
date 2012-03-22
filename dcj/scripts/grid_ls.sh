#! /bin/bash
for i in `python -c 'for i in range(1,16): print "%02d"%i,'`; do echo -e node${i}"\n\t"`ssh -fCT node${i} ls $*`; sleep 2; done
