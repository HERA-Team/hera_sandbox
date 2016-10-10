#$ -S /bin/bash
#$ -N casa_uv2ms
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
#$ -l h_vmem=4000M
CAL=psa455_v004_gc
unset DISPLAY
ARGS=`pull_args.py $*`
echo casapy --logfile grid_output/casalog_bash_uv2ms_${JOB_ID}_${TASK_ID}.log --nologger -c bash_uv2ms.py -C ${CAL} ${ARGS}
casapy --logfile grid_output/casalog_bash_uv2ms_${JOB_ID}_${TASK_ID}.log --nologger -c bash_uv2ms.py -C ${CAL} ${ARGS}
