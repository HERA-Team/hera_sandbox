#Grid script.
#Execute with: qsub -t 1:$NFILES:1 batch_srcsub.sh

DATE=`date`
EXP_NAME=srcsub
FILES=/data1/paper/arp/pgb966/lst_v002/*uvddddm
# -t 1:$TMAX:1
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
FILE=`python -c "import glob; print glob.glob('/data1/paper/arp/pgb966/lst_v002/*uvddddm')[$SGE_TASK_ID]"`
echo Working on $FILE
mdlvis.py -C pgb966_v003 -s her $FILE
