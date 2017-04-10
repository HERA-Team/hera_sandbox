#$ -S /bin/bash
#$ -j y
#$ -N rfi_surv
#$ -o grid_output/
#$ -V
#$ -cwd
ARGS=`pull_args.py $*`
rfi_sum.py $ARGS
