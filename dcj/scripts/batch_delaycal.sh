#$ -S /bin/bash
#$ -N delaycal
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
#$ -l h_vmem=4000M
CAL=psa455_v004_gc
unset DISPLAY
ARGS=`pull_args.py $*`
PREFIX=~/scripts/
echo casapy --logfile grid_output/casalog_delaycal_${JOB_ID}_${TASK_ID}.log --nologger -c ${PREFIX}/bash_delaycal.py \
-C ${CAL} -s all --cat=southern_sky_v3 --useflagversion=flagbadbls ${ARGS}

casapy --logfile grid_output/casalog_delaycal_${JOB_ID}_${TASK_ID}.log --nologger -c ${PREFIX}/bash_delaycal.py \
-C ${CAL} -s all --cat=southern_sky_v3 --useflagversion=flagbadbls ${ARGS}

