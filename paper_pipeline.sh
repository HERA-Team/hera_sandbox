#! /bin/bash
#CAL=pgb966_v003
CAL=pgb564_ver006
SUBSRCS=cyg,cas,vir,crab,Sun
#SUBSRCS=`python -c "import $CAL ; print ','.join(${CAL}.src_prms.keys())"`
FLTSRCS="Sun cyg cas vir crab"
NCHAN=256
XRFI_CH=0_239,720_1023
POL=yy
DLYWID=5
DRTWID=5
CLN=1e-3
MINUV=20

# ---------------------------------------

#FITSRCS=`python -c "print '${SRCS}'.replace(',',' ')"`

echo $SUBSRCS
for FILE in $*; do
    echo Working on $FILE
    if ! ls ${FILE}drxam &> /dev/null; then
        echo "File ${FILE}xrm doesn't exist"
        echo "Creating it..."
        mdlvis.py -C $CAL -s $SUBSRCS $FILE
        uv_addsub.py --sub $FILE ${FILE}s
        xrfi.py -m val -c $XRFI_CH -n 2 ${FILE}d
        xtalk3.py ${FILE}dr
        uv_addsub.py ${FILE}drx ${FILE}s
        rm -r ${FILE}{s,d,dr,drx}
    fi

    ARG=${FILE}drxa
    OARG=$ARG
    echo MDLVIS: $SUBSRCS
    mdlvis.py -C $CAL -s $SUBSRCS -m sub $ARG
    NARG=${ARG}s
    for SRC in $FLTSRCS ; do
        echo Isolating $SRC
        if ! ls ${NARG}*.bm_${SRC} &> /dev/null; then
            echo DDR filter: extracting $SRC residuals
            filter_src.py -p -C $CAL -e -s $SRC -d $DLYWID -r $DRTWID --clean=$CLN $NARG
            echo MDLVIS: adding $SRC back in
            mdlvis.py -C $CAL -s $SRC -m add ${NARG}.e_${SRC}
            echo Beamforming on $SRC
            beamform.py -f -C $CAL -p yy --minuv=$MINUV -s $SRC ${NARG}.e_${SRC}s
            uv_addsub.py --sub $OARG ${NARG}.e_${SRC}s
            OARG=${OARG}d
            rm -r ${NARG}.e_${SRC}
        fi
    done
    combine_freqs.py -n $NCHAN $OARG

done
