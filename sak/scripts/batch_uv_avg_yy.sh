#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=1G
#$ -j y
#$ -N uv_avg_yy
#$ -o grid_output
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Polar/bin/activate
POL='yy'
FILES=`pull_args.py $*`
for FILE in $FILES
do
echo Processing $FILE ${POL}
time python uv_avg_SK.py $FILE -p ${POL}
done
