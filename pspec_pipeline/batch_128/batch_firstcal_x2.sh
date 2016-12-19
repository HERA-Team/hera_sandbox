#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=15G
#$ -j y
#$ -N fcalx2all
#$ -o grid_output

source activate PAPER

# Hand this script a list of JDs, corresponding to the PSA-128 file structure of
# /data4/paper/2013EoR/{integer JD}/*{xx,yy}.uvcRRE

CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
JD_DIRS=`${PATH2CAPO}/scripts/pull_args.py $*`

#UBLS="64_48,64_29,64_28,64_55,64_34,64_27,64_51,64_57,64_24"
# These are the longer East-West baseline-types in the PSA-128 grid.
# Pick longer ones to hone-in on point sources in delay space (315 redundant bls).
UBLS="64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_28,64_55,64_34,64_27,64_51,64_57,64_24"


#####

for DIR in $JD_DIRS
do

OUTPATH=${DIR}

echo Working on data in ${OUTPATH}

for POL in xx yy
do

echo polarization ${POL}

FILES=${DIR}/*${POL}.uvcRRE

echo "Firstcal iteration 1"
FCAL_1_PATH=${OUTPATH}/fcal_${POL}_1_allEW
echo mkdir $FCAL_1_PATH
mkdir $FCAL_1_PATH

for FILE in $FILES
do

echo Firstcal-ing $FILE
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --outpath=$FCAL_1_PATH

done #done with single-pol firstcal-per-file loop

#1st analysis run, now that we have all firstcal files for that JD and pol
echo ${PATH2CAPO}/omni/analyze_firstcal.py --badantsout=${FCAL_1_PATH}/badants_1.txt --statfile=${FCAL_1_PATH}/stats_1.txt ${FCAL_1_PATH}/*${POL}*npz
${PATH2CAPO}/omni/analyze_firstcal.py --badantsout=${FCAL_1_PATH}/badants_1.txt --statfile=${FCAL_1_PATH}/stats_1.txt ${FCAL_1_PATH}/*${POL}*npz

#I don't understand why this is required, but it is
echo chmod 777 $FCAL_1_PATH/*
chmod 777 $FCAL_1_PATH/*
#

echo "Firstcal iteration 2"
FCAL_2_PATH=${OUTPATH}/fcal_${POL}_2_allEW
echo mkdir $FCAL_2_PATH
mkdir $FCAL_2_PATH

for FILE in $FILES #the data files have not changed name; we are still in the same POL loop
do

echo Firstcal-ing $FILE
echo ${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_1_PATH}/badants_1.txt 
EX_ANTS=`${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_1_PATH}/badants_1.txt`
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --ex_ants=${EX_ANTS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --outpath=${FCAL_2_PATH} --ex_ants=${EX_ANTS}

done #now have all firstcal files for a single pol, a second time

#Time to analyze them and get a badant list
echo ${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_1_PATH}/badants_1.txt --badantsout=${FCAL_2_PATH}/badants_2.txt --statfile=${FCAL_2_PATH}/stats_2.txt ${FCAL_2_PATH}/*${POL}*npz
${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_1_PATH}/badants_1.txt --badantsout=${FCAL_2_PATH}/badants_2.txt --statfile=${FCAL_2_PATH}/stats_2.txt ${FCAL_2_PATH}/*${POL}*npz

echo chmod 777 $FCAL_2_PATH/*
chmod 777 $FCAL_2_PATH/*

### Begin extra flagging loop

echo "Firstcal iteration 3"
FCAL_3_PATH=${OUTPATH}/fcal_${POL}_3_allEW
echo mkdir $FCAL_3_PATH
mkdir $FCAL_3_PATH

for FILE in $FILES #the data files have not changed name; we are still in the same POL loop
do

echo Firstcal-ing $FILE
echo ${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_2_PATH}/badants_2.txt 
EX_ANTS=`${PATH2CAPO}/pspec_pipeline/getExAnts.py $FILE -f ${FCAL_2_PATH}/badants_2.txt`
echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --ex_ants=${EX_ANTS}
${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} $FILE --ubls=${UBLS} --outpath=${FCAL_3_PATH} --ex_ants=${EX_ANTS}

done #now have all firstcal files for a single pol, a third time

#Time to analyze them and get a final badant list
echo ${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_2_PATH}/badants_2.txt --badantsout=${FCAL_3_PATH}/badants_3.txt --statfile=${FCAL_3_PATH}/stats_3.txt ${FCAL_3_PATH}/*${POL}*npz

${PATH2CAPO}/omni/analyze_firstcal.py --badantsin=${FCAL_2_PATH}/badants_2.txt --badantsout=${FCAL_3_PATH}/badants_3.txt --statfile=${FCAL_3_PATH}/stats_3.txt ${FCAL_3_PATH}/*${POL}*npz

echo chmod 777 $FCAL_3_PATH/*
chmod 777 $FCAL_3_PATH/*

done #done with firstcal and analysis for both pols

done #done firstcal and analysis for all JDs

echo Done-zo!

