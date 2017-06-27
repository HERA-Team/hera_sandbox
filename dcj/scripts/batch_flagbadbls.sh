#$ -S /bin/bash
#$ -N casa_flagbadbls
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
#$ -l h_vmem=8G
ARGS=`pull_args.py $*`
unset DISPLAY
echo casapy --logfile grid_output/casalog_bash_flagbadbls_${JOB_ID}_${SGE_TASK_ID}.log --nologger -c ${CASASCRIPTS}/bash_flagbadbls.py\
--cstart=300 --cstop=1700 ${ARGS}
casapy --logfile grid_output/casalog_bash_flagbadbls_${JOB_ID}_${SGE_TASK_ID}.log --nologger -c ${CASASCRIPTS}/bash_flagbadbls.py --cstart=300 --cstop=1700 ${ARGS}
