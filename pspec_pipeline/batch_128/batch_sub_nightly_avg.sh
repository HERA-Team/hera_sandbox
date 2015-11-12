#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=1G
#$ -j y
#$ -N sub_uv_xx
#$ -o grid_output
CAL=psa6942_v000
XTALK_DIR=./
#NOTE: To be run in ~jacobsda/storage/psa128/uv_avg
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
FILES=`~/scripts/pull_args.py $*`
for FILE in $FILES
do
echo Processing $FILE
time ~/scripts/sub_nightly_avg.py --xtalk_dir=${XTALK_DIR} $FILE
done

