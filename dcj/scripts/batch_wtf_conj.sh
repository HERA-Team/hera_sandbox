#$ -S /bin/bash
#$ -N psa332_C
#$ -j y
#$ -o grid_output/
ARGS=`pull_args.py $*`
~/scripts/psa331_wtf_conj.py $ARGS