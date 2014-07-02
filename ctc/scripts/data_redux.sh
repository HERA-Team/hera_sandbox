#! /bin/bash

#UV=$*
#for FILE in $UV; do
#    correct_psa898.py $FILE
#    xrfi_simple.py -a 1 --combine -t 10 -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 ${FILE}c
##    xrfi_simple.py -a all -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 ${FILE}c
#    rm -rf ${FILE}c
#done
#
#UVCR=`python -c "print ' '.join([f+'cR' for f in '${UVC}'.split()])"`

UVCR=$*
for FILE in $UVCR; do
    echo -------------------------------------
    echo Working on $FILE
    TRIPLET=`get_uv_neighbor.py $FILE`
    TRIP_LEN=`python -c "print len('${TRIPLET}'.split())"`
    echo Using files: $TRIPLET
    if [ $TRIP_LEN -lt 3 ] ; then 
        echo No adjacent files to use.  Skipping...
        exit
        continue
    fi
    TRIPLETR=`python -c "t = '${TRIPLET}'.split(); print ' '.join([t[0],t[1]+'R',t[2]])"`
    if ! ls ${FILE}R &> /dev/null; then
        echo ${FILE}R not found.  Assuming improved flags need to be generated...
        ddr_filter_coarse.py -a 1 -p xx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET || exit
        xrfi_simple.py -a 1 --combine -t 40 -n 4 ${FILE}E --to_npz=${FILE}E.npz || exit
        xrfi_simple.py -a all --combine -t 40 ${FILE} --from_npz=${FILE}E.npz || exit
    fi
    if ls ${FILE}R &> /dev/null; then
        ddr_filter_coarse.py -a all -p xx,yy,xy,yx --clean=1e-4 --maxbl=300 --nsections=10 $TRIPLETR || exit
    fi
done
