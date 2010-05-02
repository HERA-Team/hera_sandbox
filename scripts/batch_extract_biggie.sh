#$ -S /bin/bash
CAL=psa112_v003
SRCS="Sun cyg crab"
CATS=helm,misc,culgoora
RMSRCS=`python -c "print '$SRCS'.replace(' ',',')"`
ARGS=`pull_args.py $*`
DLYWID=5
DRTWID=5
CLN=1e-3
MINUV=20

for ARG in $ARGS ; do
    echo Processing $ARG
    echo MDLVIS: $RMSRCS
    mdlvis.py -C $CAL -s $RMSRCS --cat=$CATS -m sub $ARG
    #echo DDR filter: Removing Sun residuals
    #filter_src.py -p -C $CAL -s Sun -d $DLYWID -r $DRTWID --clean=$CLN ${ARG}s
    #NARG=${ARG}s.Sun
    NARG=${ARG}s
    for SRC in $SRCS ; do
        echo Isolating $SRC
        #if ! ls ${NARG}*.bm_${SRC} &> /dev/null; then
            echo DDR filter: extracting $SRC residuals
            filter_src.py -p -C $CAL -e -s $SRC --cat=$CATS -d $DLYWID -r $DRTWID --clean=$CLN $NARG
            echo MDLVIS: adding $SRC back in
            mdlvis.py -C $CAL -s $SRC --cat=$CATS -m add ${NARG}.e_${SRC}
            #echo Beamforming on $SRC
            #beamform.py -f -C $CAL -p yy --minuv=$MINUV -s $SRC ${NARG}.e_${SRC}s
        #fi
    done
done
