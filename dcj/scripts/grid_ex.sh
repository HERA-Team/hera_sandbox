#! /bin/bash
for i in `python -c 'for i in range(1,24): print "%02d"%i,'`; do echo node${i}; ssh node${i}  $*; sleep 2; done
