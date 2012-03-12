#$ -S /bin/bash
#$ -N psa332_r
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
ARGS=`pull_args.py $*`
xrfi.py -m val --blank --ch_thresh=0.6 -n 2 -c 1_300,710_800,1680_1730,1864_1871,1535_1545,1790_2045 --val_thresh=1.3e6 $ARGS
