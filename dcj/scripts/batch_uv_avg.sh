#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=1G
#$ -j y
#$ -N uv_avg
#$ -o grid_output
CAL=psa6622_v001
#NOTE: To be run in ~jacobsda/storage/psa128/uv_avg
#canopy-PAPER_Omni
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
FILES=`~/scripts/pull_args.py $*`
for FILE in $FILES
do
echo Processing $FILE
time ~/scripts/uv_avg.py $FILE
done

