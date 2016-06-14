#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=10G
#$ -j y
#$ -N pull
#$ -o grid_output

POL=xx
CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
source activate PAPER
FILES=`${PATH2CAPO}/scripts/pull_args.py $*`

UBLS="64_48,64_29,64_28,64_55,64_34,64_27,64_51,64_57,64_24" 
#pick longer ones to hone-in on point sources in delay space (315 redundant bls)

for FILE in  $FILES
do
echo Firstcal-ing $FILE

echo ${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${PATH2CAPO}/pspec_pipeline/badant4jd.txt
EX_ANTS=`${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${PATH2CAPO}/pspec_pipeline/badant4jd.txt`
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --ex_ants=${EX_ANTS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --ex_ants=${EX_ANTS}


done
