#$ -S /bin/bash
#$ -N psa332_x
#$ -V
#$ -cwd
#$ -j y
#$ -o grid_output/
#$ -l h_vmem=0.5G
ARGS=`pull_args.py $*`
xtalk3.py -o $ARGS
#xtalk3.py -r --dt=4 $ARGS
#xtalk3.py -i $ARGS