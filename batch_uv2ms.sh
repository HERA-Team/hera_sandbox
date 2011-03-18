#$ -S /bin/bash
#$ -N psa332_r
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd

CAL=psa455_v004_gc

ARGS=`pull_args.py $*`
casapy --logfile grid_output/casalog_bash_uv2ms_${JOB_ID}_${TASK_ID}.log --nologger -c bash_uv2ms.py -C ${CAL} ${ARGS}
