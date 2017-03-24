#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=5G
#$ -j y
#$ -N correct128_s2
#$ -o grid_output

source activate PAPER

# Hand this script a list of JDs, corresponding to the PSA-128 file structure of
# /data4/paper/2013EoR/{integer JD}/*.uvcRRE

CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
JD_DIRS=`${PATH2CAPO}/scripts/pull_args.py $*`
XPOLS="26,34,38,46,50,72,110"

for DIR in $JD_DIRS
do
    JD=`basename ${DIR}`
    echo Working on data in ${DIR}
    echo  ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} ${DIR}/*xx.uvcRRE
    ${PATH2CAPO}/sak/scripts/correct_psa128_pols.py -a ${XPOLS} ${DIR}/*xx.uvcRRE
done
