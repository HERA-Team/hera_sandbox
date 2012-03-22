#! /bin/bash
#$ -N Irx_cached
cd /scratch/paper/psa113/
apply_bp.py *c 
sum_integrations.py -n 10 *b 
xrfi.py -m val -c 0_33,209_255 *I 
xtalk3.py *r
