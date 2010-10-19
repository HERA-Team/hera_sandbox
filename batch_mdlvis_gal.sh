#Grid script.
#Execute with:  qsub -t 1:10:1 batch_mdlvis_gal.sh

DATE=`date`
EXP_NAME=GSM_gen
#$ -S /bin/bash
#$ -V
#$ -N GSM_gen
#$ -cwd
#$ -j y
#$ -o grid_stdout/
REGISTERED=`grep "^$JOB_ID\$" grid_stdout/log.txt`
#echo $REGISTERED
if [ -z "$REGISTERED" ]
then
        echo $JOB_ID,$DATE,$EXP_NAME >>grid_stdout/log.txt
fi
echo $SGE_TASK_ID
START=$(($SGE_TASK_ID-1))
STOP=$SGE_TASK_ID
../scripts/mdlvis_gal.py --sim -C pgb966_v003 -m all_srcs1.fits --startjd=2453887.$START --endjd=2453888.$STOP --pol=xx --sfreq=0.120 --sdf=0.000234 --nchan=256