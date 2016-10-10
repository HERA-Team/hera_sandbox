#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -N testmir
#$ -l h_vmem=500M
FILES=`~/scripts/pull_args.py $*`
source activate PAPER
#echo $FILES
python ~/scripts/test_miriad.py $FILES
