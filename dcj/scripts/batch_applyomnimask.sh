#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -j y
#$ -N flag
#$ -o grid_output
#$ -q all.q
CAL=psa6240_v003 
FILES=`~/scripts/pull_args.py $*`
source activate PAPER
python ~/scripts/apply_omniflags.py $FILES
