#$ -S /bin/bash
#$ -N psa332_r
#$ -j y
#$ -o grid_output/
ARGS=`pull_args.py $*`
~/scripts/xrfi_simple.py -c770 --dt=3 $ARGS
