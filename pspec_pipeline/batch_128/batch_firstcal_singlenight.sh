#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=10G
#$ -j y
#$ -N firstcal_single
#$ -o grid_output

source activate PAPER

# Hand this script a list of *xx* files for a single JD.

CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
FILES=`${PATH2CAPO}/scripts/pull_args.py $*`
UBLS="64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57"
FCAL_PATH=/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_4pol_xtalk/

#S2E2 ex ants - xx and yy combined badants from documented QA cuts
EX_ANTS=7,8,34,56,84,100

for FILE in ${FILES}
do
    for POL in xx yy
    do
        FP="${FILE/xx/$POL}"
        echo ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} --ubls=${UBLS} --outpath=${FCAL_PATH} --ex_ants=${EX_ANTS} $FP
        ${PATH2CAPO}/omni/firstcal.py -C ${CAL} -p ${POL} --ubls=${UBLS} --outpath=${FCAL_PATH} --ex_ants=${EX_ANTS} $FP
    done 
done
