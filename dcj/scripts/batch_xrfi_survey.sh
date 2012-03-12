#$ -S /bin/bash 
#$ -j y
#$ -N rfi_surv
#$ -cwd
#$ -V
ARGS=`pull_args.py -w $*`
echo xrfi_filemaker.py -n 2.5 -m val -o $ARGS
xrfi_filemaker.py -n 2.5 -m val -o $ARGS
