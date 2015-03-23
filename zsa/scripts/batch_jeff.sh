#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G

ARGS=`pull_args.py $*`
for dir in $ARGS; 
