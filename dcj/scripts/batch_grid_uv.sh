#$ -S /bin/bash
#$ -N grid_uv
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
#$ -l h_vmem=3G
unset DISPLAY
ARGS=`pull_args.py $*`

UVSIZE=200
UVRES=4
CHANS=300_450 
echo casapy --logfile grid_output/casalog_selfcal_${JOB_ID}_${SGE_TASK_ID}.log\
 --nologger -c ~/scripts/bash_grid_uv.py \
--uvsize=${UVSIZE} --uvres=${UVRES} --chan=${CHANS} \
$ARGS

casapy --logfile grid_output/casalog_selfcal_${JOB_ID}_${SGE_TASK_ID}.log\
 --nologger -c ~/scripts/bash_grid_uv.py \
--uvsize=${UVSIZE} --uvres=${UVRES} --chan=${CHANS} \
$ARGS

