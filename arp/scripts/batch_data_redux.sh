#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
ARGS=`pull_args.py $*`
echo data_redux.sh $ARGS
STARTPATH=`pwd`
SCRATCH=/scratch/paper

for FILE in $ARGS; do
    echo --------------------
    echo --- Working on $FILE ---
    # Only process this file if the end product is missing
    if ls ${FILE}cRRE &> /dev/null; then
        echo ${FILE}cRRE exists.  Skipping...
        continue
    fi
    FILEDIR=`python -c "import os; print os.path.dirname('$FILE')"`
    FILEBASE=`python -c "import os; print os.path.basename('$FILE')"`
    # Copy files to a local directory for fast access
    for TFILE in `get_uv_neighbor.py $FILE`; do
        TFILEBASE=`python -c "import os; print os.path.basename('$TFILE')"`
        # If the uvcR file is here already, don't bother
        if ls ${SCRATCH}/${TFILEBASE}cR &> /dev/null; then
            echo Using ${SCRATCH}/${TFILEBASE}cR
            continue
        fi
        echo Generating ${SCRATCH}/${TFILEBASE}cR
        echo cp -r $TFILE ${SCRATCH}
        cp -r $TFILE ${SCRATCH}
        correct_psa898.py ${SCRATCH}/$TFILEBASE
        #correct_psa746_v002.py -t . ${SCRATCH}/$TFILEBASE
        rm -rf ${SCRATCH}/${TFILEBASE}
        xrfi_simple.py -a 1 --combine -t 20 -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6  ${SCRATCH}/${TFILEBASE}c
        rm -rf ${SCRATCH}/${TFILEBASE}c
        echo
    done
    echo cd ${SCRATCH}
    cd ${SCRATCH}
    echo data_redux.sh ${FILEBASE}cR
    data_redux.sh ${FILEBASE}cR
    echo cd $STARTPATH
    cd $STARTPATH
    echo cp -r ${SCRATCH}/${FILEBASE}cRR[DEF] $FILEDIR
    cp -r ${SCRATCH}/${FILEBASE}cRR[DEF] $FILEDIR
    echo cp ${SCRATCH}/${FILEBASE}cRE.npz $FILEDIR
    cp ${SCRATCH}/${FILEBASE}cRE.npz $FILEDIR
    rm -rf ${SCRATCH}/${TFILE}cRR
done

# Final clean-up
for FILE in $ARGS; do
    for TFILE in `get_uv_neighbor.py $FILE`; do
        FILEBASE=`python -c "import os; print os.path.basename('$TFILE')"`
        echo rm -rf ${SCRATCH}/${FILEBASE}*
        rm -rf ${SCRATCH}/${FILEBASE}*
    done
done

echo DONE
