#! /bin/bash
CAL=psa112_v005
SUBSRCS=Sun,cyg,crab,hyd,pic,for,cen,J1615-605,J1935-461,J2154-692
#SUBSRCS=`python -c "import $CAL ; print ','.join(${CAL}.src_prms.keys())"`
#FLTSRCS="Sun cyg crab cen hyd pic for"
FLTSRCS="Sun cyg crab"
CATS=helm,misc,culgoora
NCHAN=256
XRFI_CH=0_49,200_255
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
    #if ! ls ${FILE}drxam &> /dev/null; then
    #    echo "File ${FILE}xrm doesn't exist"
    #    echo "Creating it..."
    #    #mdlvis.py -C $CAL -s $SUBSRCS $FILE
    #    #uv_addsub.py --sub $FILE ${FILE}s
    #    #xrfi.py -m val -c $XRFI_CH -n 2.5 ${FILE}d
    #    xrfi.py -m val -c $XRFI_CH -n 2.5 ${FILE}
    #    #xtalk3.py ${FILE}dr
    #    #uv_addsub.py ${FILE}drx ${FILE}s
    #    #rm -r ${FILE}{s,d,dr,drx}
    #fi

    ##ARG=${FILE}drxa
    ARG=${FILE}
    #OARG=$ARG
    #echo MDLVIS: $SUBSRCS
    #mdlvis.py -C $CAL -s $SUBSRCS --cat=$CATS -m sub $ARG
    #NARG=${ARG}s
    #uv_addsub.py --sub ${FILE} ${FILE}s
    #xtalk3.py ${FILE}d
    NARG=${ARG}dx
    OARG=${NARG}
    for SRC in $FLTSRCS ; do
        echo Isolating $SRC
        #if ! ls ${NARG}*.bm_${SRC} &> /dev/null; then
            echo DDR filter: extracting $SRC residuals
            filter_src.py -p -C $CAL -e -s $SRC --cat=$CATS -d $DLYWID -r $DRTWID --clean=$CLN $NARG
            #echo MDLVIS: adding $SRC back in
            #mdlvis.py -C $CAL -s $SRC -m add ${NARG}.e_${SRC}
            #echo Beamforming on $SRC
            #beamform.py -f -C $CAL -p yy --minuv=$MINUV -s $SRC ${NARG}.e_${SRC}s
            #uv_addsub.py --sub $OARG ${NARG}.e_${SRC}s
            uv_addsub.py --sub $OARG ${NARG}.e_${SRC}
            OARG=${OARG}d
            #rm -r ${NARG}.e_${SRC}
        #fi
    done
    #combine_freqs.py -n $NCHAN $OARG

done
