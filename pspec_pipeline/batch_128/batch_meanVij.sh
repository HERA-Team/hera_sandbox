#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=8G
#$ -j y
#$ -N meanVij
#$ -o grid_output

source activate PAPER

# Hand this script a list of JDs, corresponding to the PSA-128 file structure of
# /data4/paper/2013EoR/{integer JD}/*.uvcRRE

CAL=psa6622_v003
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
JD_DIRS=`${PATH2CAPO}/scripts/pull_args.py $*`
FX2ANTS="97,10,91,58,72,107,105,64,42,31,43,33,15,22,47,2"
OUTPATH=/home/saulkohn/flags_epoch2

for DIR in $JD_DIRS
do
    JD=`basename ${DIR}`
    echo Working on data in ${DIR}
    
    for POL in xx yy
    do
        if (( ${JD} < '2456677'))
            then
                echo ${PATH2CAPO}/sak/scripts/meanVij.py -C ${CAL} -p ${POL} --ba=${FX2ANTS} ${DIR}/*${POL}.uvcRREc --outpath=${OUTPATH}/${JD}.${POL}.txt --skiplast
                ${PATH2CAPO}/sak/scripts/meanVij.py -C ${CAL} -p ${POL} --ba=${FX2ANTS} ${DIR}/*${POL}.uvcRREc --outpath=${OUTPATH}/${JD}.${POL}.txt --skiplast
            else
                echo ${PATH2CAPO}/sak/scripts/meanVij.py -C ${CAL} -p ${POL} ${DIR}/*${POL}.uvcRRE --outpath=${OUTPATH}/${JD}.${POL}.txt --skiplast
                ${PATH2CAPO}/sak/scripts/meanVij.py -C ${CAL} -p ${POL} ${DIR}/*${POL}.uvcRRE --outpath=${OUTPATH}/${JD}.${POL}.txt --skiplast
        fi
    done
done
