#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=30G
#$ -N OMNI_RUN_4POL
source activate PAPER
FILES=`pull_args.py $*` #hand this script a single pol's-worth of files
PATH2CAPO='/home/saulkohn/ReposForCanopy/capo/'
CAL='psa6622_v003'
POL='xx,xy,yx,yy'
#OMNIPATH='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_4pol_Vmin/'
OMNIPATH='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_4pol_xtalk/'
FCALFILES='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_4pol_xtalk/fcals/zen.2456680.17382.xx.uvcRREc.median.fc.npz,/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_4pol_xtalk/fcals/zen.2456680.17382.yy.uvcRREc.median.fc.npz'
BADANTS='7,8,34,56,84,100' #these are the union of S1E2 xx and yy badants from pre-omnical QA

for FILE in ${FILES}
    do
    if (( ${FILE:33:7} > 2456678 )) #S1E2
        then #forgive me for hardcoding for folios filepaths there... correct way would be a regex
            echo ${PATH2CAPO}/omni/omni_run.py -p ${POL} -C ${CAL} --omnipath=${OMNIPATH} --fc2=${FCALFILES} --ba=${BADANTS} ${FILE}
            ${PATH2CAPO}/omni/omni_run.py -p ${POL} -C ${CAL} --omnipath=${OMNIPATH} --fc2=${FCALFILES} --ba=${BADANTS} ${FILE}
    else
        echo ${FILE} in S1E1. Skipping...
    fi
done


