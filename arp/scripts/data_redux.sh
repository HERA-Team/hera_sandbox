#! /bin/bash

UV=$*
for FILE in $UV; do
    correct_psa898.py $FILE
    xrfi_simple.py -a all -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 ${FILE}c
done

UVCR=`python -c "print ' '.join([f+'cR' for f in '${UVC}'.split()])"`
#UVCR=$*
# XXX is the "pol" argument in ddr_filter_coarse necessary?
for FILE in $UVCR; do
    echo -------------------------------------
    echo Working on $FILE
    TRIPLET=`get_uv_neighbor.py $FILE`
    TRIP_LEN=`python -c "print len('${TRIPLET}'.split())"`
    if [ $TRIP_LEN -lt 3 ] ; then 
        echo No adjacent files to use.  Skipping...
        continue
    fi
    if ! ls ${FILE}/flags_bk &> /dev/null; then
        echo ${FILE}/flags_bk not found.  Assuming improved flags need to be generated...
        ddr_filter_coarse.py -p xx -a all --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET
        xrfi_simple.py -a all -n 4 ${FILE}E
        # Check that we got to here successfully before continuing
        if [ $? -eq 0 ] ; then
            echo Copying ${FILE}ER/flags to ${FILE}
            cp ${FILE}/flags ${FILE}/flags_bk
            cp ${FILE}ER/flags ${FILE}/flags
            rm -rf ${FILE}E
        fi
        #rm -rf ${FILE}ER
    fi
    # Make sure flags_bk exists (a sign the above completed successfully) before continuing
    if ls ${FILE}/flags_bk &> /dev/null; then
        ddr_filter_coarse.py -p xx -a all --clean=1e-4 --maxbl=300 $TRIPLET
    fi
done
