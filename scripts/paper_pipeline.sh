#! /bin/bash
SRCS=cyg,cas,vir,crab,Sun
LOC=pgb564
NCHAN=256
XRFI_CH=0_239,720_1023
CH=60_180
POL=yy
XTOL=1e-3
MAXITER=200
DLYWID=3
FNGWID=3
CLN=1e-3

# ---------------------------------------

FITSRCS=`echo $SRCS | python -c "import sys; print sys.stdin.read().replace(',',' ')"`

for FILE in $*; do
    echo Working on $FILE
    JD=`echo $FILE | python -c "import sys; print '.'.join(sys.stdin.read().split('.')[1:3])"`
    UVCB=zen.${JD}.uvcb
    if ! ls ${UVCB}xrm &> /dev/null; then
        echo "File ${UVCB}xrm doesn't exist"
        echo "Creating it..."
        mdlvis.py -l $LOC -s $SRCS $UVCB
        xrfi.py -c $XRFI_CH ${UVCB}s
        xtalk3.py -o ${UVCB}sr
        xtalk3.py -i ${UVCB}s
        xrfi.py -c $XRFI_CH -o ${UVCB}sx
        xtalk3.py -i ${UVCB}
        xrfi.py -c $XRFI_CH -i ${UVCB}x
        combine_freqs.py -n $NCHAN ${UVCB}xr
        rm -r ${UVCB}{s,sr,sx,x}
    fi

    for SRC in $FITSRCS; do
        if ! ls ${JD}.${SRC} &> /dev/null; then
            echo "Ionospheric calibration file ${JD}.${SRC} is missing"
            echo "Creating it..."
            if ! ls ${JD}.${SRC}pos.fit &> /dev/null; then
                echo "Creating ${JD}.${SRC}pos.fit"
                fitmdl.py -p $POL -c $CH -l $LOC -s $SRCS --snap --sprms=${SRC}=dra,${SRC}=ddec -q --maxiter=$MAXITER --xtol=$XTOL ${UVCB}xrm | tee ${JD}.${SRC}pos.fit
            fi
            if ! ls ${JD}.${SRC}amp.fit &> /dev/null; then
                echo "Creating ${JD}.${SRC}amp.fit"
                cat ${JD}.${SRC}pos.fit | scripts/src_filter.py -l $LOC -s ${SRC} > ${JD}.${SRC}
                fitmdl.py -p $POL -c $CH -l $LOC -s $SRCS --snap --sprms=${SRC}=str,${SRC}=index -q --maxiter=$MAXITER --xtol=$XTOL ${UVCB}xrm | tee ${JD}.${SRC}amp.fit
            fi
            cat ${JD}.${SRC}pos.fit ${JD}.${SRC}amp.fit | scripts/src_filter.py -l $LOC -s ${SRC} > ${JD}.${SRC}
        fi
    done

    if ! ls ${UVCB}xrms &> /dev/null; then
        echo "${UVCB}xrms doesn't exist"
        echo "Creating it using current calibration..."
        mdlvis.py -l $LOC -s $SRCS ${UVCB}xrm
    fi

    #if ! ls ${UVCB}xrms.Sun &> /dev/null; then
    #    echo "${UVCB}xrms.Sun doesn't exist"
    #    echo "Creating it using current calibration..."
    #    filter_src.py -l $LOC -s Sun -d $DLYWID -f $FNGWID --clean=$CLN ${UVCB}xrms
    #fi
done
