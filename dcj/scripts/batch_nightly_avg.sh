#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -j y
#$ -N nightly_xtalk
#$ -o grid_output
CAL=psa6622_v001
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
JDS=`~/scripts/list_jds.py $*`
MYJDS=`~/scripts/pull_args.py ${JDS}`
echo working on JDS: $MYJDS
for JD in $MYJDS
do
    MYFILES=`python -c "print ' '.join([l for l in '''$*'''.split() if l.find('${JD}')>0])"`
    echo Working on $MYFILES
    ~/scripts/nightly_avg.py $MYFILES
done

