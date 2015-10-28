#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=1G
#$ -j y
#$ -N sub_4pol_uv
#$ -o grid_output
CAL=psa6240_v003
XTALK_DIR=night_xtalk_pkls/
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Polar/bin/activate
FILES=`pull_args.py $*`
for FILE in $FILES
do
echo Processing $FILE
time python sub_nightly_avg_SK.py --xtalk_dir=${XTALK_DIR} $FILE
done

