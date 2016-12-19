#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=5G
#$ -j y
#$ -N correct128
#$ -o grid_output

source activate PAPER

# Hand this script a list of JDs, corresponding to the PSA-128 file structure of
# /data4/paper/2013EoR/{integer JD}/*.uvcRRE

CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
JD_DIRS=`${PATH2CAPO}/scripts/pull_args.py $*`
XPOLS="26,34,38,46,50,72"
FX2ANTS="97,10,91,58,72,107,105,64,42,31,43,33,15,22,47,2"

for DIR in $JD_DIRS
do
    JD=`basename ${DIR}`
    echo Working on data in ${DIR}
    
    if (( ${JD} < '2456677'))
        then
            echo ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} -b ${FX2ANTS} ${DIR}/*xx.uvcRRE
            ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} -b ${FX2ANTS} ${DIR}/*xx.uvcRRE
        else
            echo  ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} ${DIR}/*xx.uvcRRE
            ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} ${DIR}/*xx.uvcRRE
    fi
done
