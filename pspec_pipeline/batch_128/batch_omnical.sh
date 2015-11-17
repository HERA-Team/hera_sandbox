#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=8G
#$ -j y
#$ -N omnical
#$ -o grid_output
CAL=psa6622_v001
#NOTE: To be run in ~jacobsda/storage/psa128
#NOTE: Use full path.
#. /usr/global/paper/CanopyVirtualEnvs/shredddercanopyrc.sh
#canopy-PAPER_Omni
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate OMNI
FILES=`~/scripts/pull_args.py $*`
for FILE in $FILES
do
FILEBASE=`python -c "import os; print os.path.basename('${FILE}')"`
echo Processing $FILE
echo omnical_PSA128.py -C ${CAL} -pyy -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3 -i \
redundantinfo_first_cal_epoch3yy.bin -r calpar_first_cal_epoch3yy.p -u -d $FILEBASE -t O -s $FILE 
time omnical_PSA128.py -C ${CAL} -pyy -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/ -i \
redundantinfo_first_cal_epoch3yy.bin -r calpar_first_cal_epoch3yy.p -u -d $FILEBASE -t O -s $FILE 

done

