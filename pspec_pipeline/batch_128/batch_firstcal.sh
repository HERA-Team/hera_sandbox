#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=10G
#$ -j y
#$ -N firstcal
#$ -o grid_output

source activate PAPER

## XXX TODO: these need to be wrapped into a .cfg file
POL=xx
CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
FILES=`${PATH2CAPO}/scripts/pull_args.py $*`
UBLS="64_48,64_29,64_28,64_55,64_34,64_27,64_51,64_57,64_24" 
#pick longer ones to hone-in on point sources in delay space (315 redundant bls)
OUTPATH="/data4/paper/2013EoR/2456678/"

#####

echo "Firstcal iteration 1"
FCAL_1_PATH=${OUTPATH}/fcal_1
echo mkdir $FCAL_1_PATH
mkdir $FCAL_1_PATH

for FILE in  $FILES
do
echo Firstcal-ing $FILE
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --outpath=$FCAL_1_PATH
done

echo ${PATH2CAPO}/omni/analyze_firstcal.py --badantsout=${FCAL_1_PATH}/badants_1.txt --statfile=${FCAL_1_PATH}/stats_1.txt ${FCAL_1_PATH}/*${POL}*npz
${PATH2CAPO}/omni/analyze_firstcal.py --badantsout=${FCAL_1_PATH}/badants_1.txt --statfile=${FCAL_1_PATH}/stats_1.txt ${FCAL_1_PATH}/*${POL}*npz

#####
echo chmod 777 $FCAL_1_PATH/* #WTF is this an issue
chmod 777 $FCAL_1_PATH/*
#####

echo "Firstcal iteration 2"
FCAL_2_PATH=${OUTPATH}/fcal_2
echo mkdir $FCAL_2_PATH
mkdir $FCAL_2_PATH

for FILE in  $FILES
do
echo Firstcal-ing $FILE
echo ${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_1_PATH}/badants_1.txt 
#XXX this only actually needs to be run once in the current PAPER data storage structure... could be moved out of loop
EX_ANTS=`${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_1_PATH}/badants_1.txt`
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --ex_ants=${EX_ANTS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --outpath=${FCAL_2_PATH} --ex_ants=${EX_ANTS}
done

echo ${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_1_PATH}/badants_1.txt --badantsout=${FCAL_2_PATH}/badants_2.txt --statfile=${FCAL_2_PATH}/stats_2.txt ${FCAL_2_PATH}/*${POL}*npz
${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_1_PATH}/badants_1.txt --badantsout=${FCAL_2_PATH}/badants_2.txt --statfile=${FCAL_2_PATH}/stats_2.txt ${FCAL_2_PATH}/*${POL}*npz

#####
echo chmod 777 $FCAL_2_PATH/* #WTF is this an issue
chmod 777 $FCAL_2_PATH/*
#####
